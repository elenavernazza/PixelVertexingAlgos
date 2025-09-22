import numpy as np
import ROOT
import argparse

#######################################
# Parameters to be otpimized
#######################################
ptmin = 1
ptmax = 75
minT = 2
eps = 0.07
errmax = 0.01
chi2max = 9.0
#######################################

def cluster_tracks_by_density(
    zt: np.ndarray,
    ezt2: np.ndarray,
    minT: int,
    eps: float,
    errmax: float,
    chi2max: float,
    verbose: bool = False
):
    """
    Translate of the Alpaka `clusterTracksByDensity` algorithm to Python/NumPy.

    Parameters
    ----------
    zt : array_like (nt,)
        z positions of tracks (float).
    ezt2 : array_like (nt,)
        squared errors on z (float). In C++ this is 'ezt2'.
    minT : int
        minimum number of neighbours to be a seed (density threshold).
    eps : float
        maximum absolute distance to consider two tracks neighbours.
        (C++ assumed eps <= 0.1 because of binning choice.)
    errmax : float
        maximum allowed error (not squared) to consider a track as candidate seed.
    chi2max : float
        maximum allowed normalized distance: (dist^2) <= chi2max*(ezt2_i + ezt2_j).
    verbose : bool
        print debug info.

    Returns
    -------
    iv : np.ndarray (nt,)
        cluster id per track (integers in 0..n_clusters-1). Noise will not be given a special id
        (in original code noise tracks are excluded because they get no cluster).
    nclusters : int
        number of clusters found.
    info : dict
        additional diagnostic arrays: 'nn' (neighbour counts), 'iv_intermediate' (before final relabel),
        'seed_mask', 'raw_iv_before_shift'.
    """

    nt = len(zt)
    assert ezt2.shape[0] == nt

    # Working arrays analogous to C++ code:
    # izt: bin index stored as uint8 (C++ did iz = int(zt*10), clamp to int8, then izt = iz - INT8_MIN)
    # we replicate that mapping, but in python we'll keep ints.
    iz = (zt * 10.0).astype(int)  # this is valid if eps <= 0.1 as in the original
    iz = np.clip(iz, -128, 127)   # clamp to int8
    izt = iz - (-128)             # shift to range [0,255] like uint8 in C++ (i.e. iz - INT8_MIN)
    izt = izt.astype(np.int16)    # small integer type for indexing convenience

    if verbose:
        print(f"nt={nt}, eps={eps}, errmax={errmax}, chi2max={chi2max}")

    # Build histogram-like mapping: bin -> list of indices that fall into it.
    # In C++ they used a specialized histogram container; here we use a dict of lists.
    bins: Dict[int, list] = {}
    for i in range(nt):
        b = int(izt[i])
        bins.setdefault(b, []).append(i)

    # Initialize iv (assignment pointer), nn (neighbor count)
    iv = np.arange(nt, dtype=int)   # initially iv[i] = i
    nn = np.zeros(nt, dtype=int)    # neighbor counts (like trkdata.ndof() being used as nn in C++)

    errmax2 = errmax * errmax

    # Count neighbours
    for i in range(nt):
        if ezt2[i] > errmax2:
            # skip tracks with too-large error
            continue
        b = int(izt[i])
        # check neighboring bins b-1, b, b+1 (C++ forEachInBins(..., 1, ...))
        for bb in (b-1, b, b+1):
            if bb not in bins:
                continue
            for j in bins[bb]:
                if j == i:
                    continue
                dist = abs(zt[i] - zt[j])
                if dist > eps:
                    continue
                # chi2-like check
                if dist*dist > chi2max * (ezt2[i] + ezt2[j]):
                    continue
                nn[i] += 1

    # Find closest "above me" (as in C++). The logic:
    # for each i set mdist = eps, loop candidate j in neighboring bins:
    #   skip j if nn[j] < nn[i]
    #   if nn[j] == nn[i] and zt[j] >= zt[i] -> skip (natural ordering)
    #   if dist>mdist skip
    #   if chi2 condition fails skip
    #   else mdist=dist; iv[i]=j
    for i in range(nt):
        mdist = eps
        b = int(izt[i])
        for bb in (b-1, b, b+1):
            if bb not in bins:
                continue
            for j in bins[bb]:
                if j == i:
                    continue
                if nn[j] < nn[i]:
                    continue
                if nn[j] == nn[i] and zt[j] >= zt[i]:
                    continue
                dist = abs(zt[i] - zt[j])
                if dist > mdist:
                    continue
                if dist*dist > chi2max * (ezt2[i] + ezt2[j]):
                    continue
                mdist = dist
                iv[i] = j

    # Consolidate graph: follow pointers until fixed point (root)
    # Equivalent to: for each i: m = iv[i]; while (m != iv[m]) m = iv[m]; iv[i] = m;
    for i in range(nt):
        m = iv[i]
        # path-following (no compression at first)
        while m != iv[m]:
            m = iv[m]
        iv[i] = m

    # Find seeds (tracks where iv[i] == i). If nn[i] >= minT => seed, else noise.
    foundClusters = 0
    # Use same marking as C++: seeds get iv[i] = -(old+1), noise gets -9998
    for i in range(nt):
        if iv[i] == i:
            if nn[i] >= minT:
                old = foundClusters
                iv[i] = -(old + 1)   # negative id for cluster seeds
                foundClusters += 1
            else:
                iv[i] = -9998        # noise

    # Propagate the negative id to all tracks in the cluster:
    # for i: if iv[i] >= 0 then iv[i] = iv[ iv[i] ]
    for i in range(nt):
        if iv[i] >= 0:
            iv[i] = iv[iv[i]]

    # Adjust cluster id to be positive starting from 0: iv[i] = -iv[i] - 1
    # But be careful: noise entries will be -9998 and will become positive large numbers:
    # original C++ treats noise specially (they ended up with big value). We'll mark noise as -1 afterwards.
    for i in range(nt):
        if iv[i] < 0:
            if iv[i] == -9998:
                iv[i] = -1     # mark noise as -1
            else:
                iv[i] = -iv[i] - 1

    # at this point iv in [0..foundClusters-1] for clustered tracks, -1 for noise
    nclusters = foundClusters

    if verbose:
        print(f"foundClusters (seed count) = {foundClusters}")
        # optionally print cluster sizes
        if nclusters > 0:
            uniq, counts = np.unique(iv[iv >= 0], return_counts=True)
            print("cluster sizes:", dict(zip(uniq.tolist(), counts.tolist())))
        nnoise = np.sum(iv == -1)
        print(f"noise tracks: {nnoise}")

    info = {
        "nn": nn,
        "nvIntermediate": foundClusters,
        "raw_iv_after_consolidation": iv.copy(),  # final iv (with -1 for noise)
    }

    return iv, nclusters, info

def fit_vertices(
    zt: np.ndarray,
    ezt2: np.ndarray,
    iv: np.ndarray,
    nclusters: int,
    chi2Max: float,
    verbose: bool = False,
):
    """
    Python translation of the Alpaka `fitVertices` kernel.
    
    Parameters
    ----------
    zt : np.ndarray (nt,)
        z positions of tracks
    ezt2 : np.ndarray (nt,)
        squared z errors of tracks
    iv : np.ndarray (nt,)
        cluster id assignment for tracks (from cluster_tracks_by_density),
        -1 = noise
    nclusters : int
        number of clusters found
    chi2Max : float
        maximum per-track chi2 contribution to accept the track
    verbose : bool
        print debug info
    
    Returns
    -------
    zv : np.ndarray (nclusters,)
        fitted vertex z positions
    wv : np.ndarray (nclusters,)
        fitted vertex weights
    chi2 : np.ndarray (nclusters,)
        per-vertex chi2
    iv_out : np.ndarray (nt,)
        updated cluster assignment (tracks rejected by chi2 get iv=9999)
    info : dict
        debug info including ndof per vertex
    """
    nt = len(zt)
    assert len(ezt2) == nt
    assert len(iv) == nt

    # outputs
    zv = np.zeros(nclusters, dtype=float)
    wv = np.zeros(nclusters, dtype=float)
    chi2 = np.zeros(nclusters, dtype=float)
    pt2 = np.zeros(nclusters, dtype=float)
    nn = np.full(nclusters, -1, dtype=int)  # degrees of freedom per vertex

    iv_out = iv.copy()

    # --- compute cluster location (weighted average) ---
    for i in range(nt):
        if iv_out[i] < 0:   # noise tracks
            continue
        if iv_out[i] >= nclusters:
            continue
        w = 1.0 / ezt2[i]
        zv[iv_out[i]] += zt[i] * w
        wv[iv_out[i]] += w

    for k in range(nclusters):
        if wv[k] <= 0:
            if verbose:
                print(f"WARNING: vertex {k} has no weight!")
            continue
        zv[k] /= wv[k]
        nn[k] = -1  # reset ndof

    # --- compute chi2 and outlier rejection ---
    for i in range(nt):
        if iv_out[i] < 0 or iv_out[i] >= nclusters:
            continue
        c2 = (zv[iv_out[i]] - zt[i]) ** 2 / ezt2[i]
        if c2 > chi2Max:
            iv_out[i] = 9999  # mark as outlier
            continue
        chi2[iv_out[i]] += c2
        if nn[iv_out[i]] < 0:
            nn[iv_out[i]] = 0
        nn[iv_out[i]] += 1

    # --- adjust weights ---
    for k in range(nclusters):
        if nn[k] > 0 and chi2[k] > 0:
            wv[k] *= float(nn[k]) / chi2[k]

    info = {
        "ndof": nn,
        "nclusters": nclusters,
        "nnoise": np.sum(iv < 0),
        "noutliers": np.sum(iv_out == 9999),
    }

    if verbose:
        print(f"fit_vertices: {nclusters} clusters")
        print("zv:", zv)
        print("wv:", wv)
        print("chi2:", chi2)
        print("ndof:", nn)
        print("noise:", info["nnoise"], "outliers:", info["noutliers"])

    return zv, wv, chi2, nn, iv_out, info


def map_vertices(my_z, my_z_err, my_ntracks, ref_z, ref_z_err, ref_ntracks, max_delta=0.05):
    """
    Map fitted vertices (my_z) to reference vertices (ref_z) by nearest z.

    Parameters
    ----------
    my_z : array-like
        Fitted vertex z positions
    my_z_err : array-like
        Errors on fitted z
    my_ntracks : array-like
        Number of tracks in each fitted vertex
    ref_z : array-like
        Reference vertex z positions
    ref_z_err : array-like
        Errors on reference z
    ref_ntracks : array-like
        Reference vertex track counts
    max_delta : float
        Maximum |z_my - z_ref| allowed for a match (in cm)

    Returns
    -------
    matches : list of dict
        Each entry has keys:
          "my_z", "my_z_err", "my_ntracks",
          "ref_z", "ref_z_err", "ref_ntracks",
          "delta_z", "delta_tracks"
    unmatched_my : list of dict
        Fitted vertices that could not be matched
    unmatched_ref : list of dict
        Reference vertices that were not matched
    """
    my_z = np.asarray(my_z)
    my_z_err = np.asarray(my_z_err)
    my_ntracks = np.asarray(my_ntracks)
    ref_z = np.asarray(ref_z)
    ref_z_err = np.asarray(ref_z_err)
    ref_ntracks = np.asarray(ref_ntracks)

    matches = []
    unmatched_my = []
    unmatched_ref = []

    used_ref = set()

    for i, (mz, mz_err, mnt) in enumerate(zip(my_z, my_z_err, my_ntracks)):
        deltas = np.abs(ref_z - mz)
        j = np.argmin(deltas)
        if deltas[j] < max_delta:
            matches.append({
                "my_z": mz, "my_z_err": mz_err, "my_ntracks": mnt,
                "ref_z": ref_z[j], "ref_z_err": ref_z_err[j], "ref_ntracks": ref_ntracks[j],
                "delta_z": deltas[j], "delta_tracks": mnt - ref_ntracks[j]
            })
            used_ref.add(j)
        else:
            unmatched_my.append({
                "my_z": mz, "my_z_err": mz_err, "my_ntracks": mnt
            })

    for j, (rz, rz_err, rnt) in enumerate(zip(ref_z, ref_z_err, ref_ntracks)):
        if j not in used_ref:
            unmatched_ref.append({
                "ref_z": rz, "ref_z_err": rz_err, "ref_ntracks": rnt
            })

    return matches, unmatched_my, unmatched_ref

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Make HLT Jet validation plots.')
    parser.add_argument('-f', '--file', type=str, required=True, default='step2_L1P2GT_HLT_DQM_VALIDATION_NANO_inNANOAODSIM.root', help='Paths to the DQM ROOT file.')
    args = parser.parse_args()

    # Read input file and create RDataFrame
    files = [args.file]
    dataframe_files = ROOT.vector(str)()
    for f in files:
        file = ROOT.TFile.Open(f)
        if file.Get("Events"):
            dataframe_files.push_back(f)
    df = ROOT.RDataFrame("Events", dataframe_files)

    # Read relevant branches
    entries = df.AsNumpy(columns=["hltPixelTrack_pt", "hltPixelTrack_dZ", "hltPixelTrack_dZError", "hltPhase2PixelVertices_z", "hltPhase2PixelVertices_dz", "hltPhase2PixelVertices_tracksSize"])

    # Read tracks and produce vertices with clustering algorithm
    for i, (pts, zs, dzs) in enumerate(zip(entries["hltPixelTrack_pt"], 
                                           entries["hltPixelTrack_dZ"], 
                                           entries["hltPixelTrack_dZError"])):
        
        if i == 0:

            # remove tracks with pt outside [ptmin, ptmax]
            pts = np.array(pts, dtype=float)
            zs  = np.array(zs, dtype=float)
            dzs = np.array(dzs, dtype=float)
            mask = (pts >= ptmin) & (pts <= ptmax)
            zs = zs[mask]
            dzs = dzs[mask]

            # cluster vertices
            iv, nclusters, info = cluster_tracks_by_density(zs, np.square(dzs), minT, eps, errmax, chi2max, verbose=False)
            zv, wv, chi2, nn, iv_out, info = fit_vertices(zs, np.square(dzs), iv, nclusters, chi2Max=50.0, verbose=False)
            print("My clustering results:")
            for i in range(len(zv)): 
                print(f"Vertex {i}: z[{i}] = {zv[i]}, dz[{i}] = {1./np.sqrt(wv[i])}, ntrks[{i}] = {nn[i]}")
            my_z = zv
            my_z_err = 1./np.sqrt(wv)
            my_ntracks = nn


    # Compare with reference vertices from NanoAOD
    for i, (zs, dzs, trks) in enumerate(zip(entries["hltPhase2PixelVertices_z"], 
                                            entries["hltPhase2PixelVertices_dz"],
                                            entries["hltPhase2PixelVertices_tracksSize"])):
        if i == 0:
            print("Reference vertices from NanoAOD:")
            for i in range(len(zs)):
                print(f"Vertex {i}: z[{i}] = {zs[i]}, dz[{i}] = {dzs[i]}, ntrks[{i}] = {trks[i]}")
            ref_z = zs
            ref_z_err = dzs
            ref_ntracks = trks

    matches, unmatched_my, unmatched_ref = map_vertices(
        my_z, my_z_err, my_ntracks,
        ref_z, ref_z_err, ref_ntracks,
        max_delta=0.01
    )

    for m in matches:
        print(f"My z={m['my_z']:.3f}±{m['my_z_err']:.4f} ({m['my_ntracks']} trks) "
            f"↔ Ref z={m['ref_z']:.3f}±{m['ref_z_err']:.4f} ({m['ref_ntracks']} trks), "
            f"Δz={m['delta_z']:.3f}, Δtracks={m['delta_tracks']}")
    for um in unmatched_my:
        print(f"My unmatched vertex: z={um['my_z']:.3f}±{um['my_z_err']:.4f} ({um['my_ntracks']} trks)")
    for ur in unmatched_ref:
        print(f"Ref unmatched vertex: z={ur['ref_z']:.3f}±{ur['ref_z_err']:.4f} ({ur['ref_ntracks']} trks)")
    print(f"Total matches: {len(matches)}, unmatched my: {len(unmatched_my)}, unmatched ref: {len(unmatched_ref)}")
