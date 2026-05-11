#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import MultiSuSiE
import os

# ----------------------------
# CONFIG
# ----------------------------

POPULATIONS = ["AFR", "AMR", "EAS", "EUR", "SAS", "UGR"]

FILE_PATHS = {
    "AFR": "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/ldlc/ukb_afr_ldlc_common_snps_formatted.txt",
    "AMR": "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/ldlc/mcps_ldlc_common_snps_formatted.txt",
    "EAS": "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/ldlc/ckb_ldlc_common_snps_formatted.txt",
    "EUR": "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/ldlc/ukb_eur_ldlc_common_snps_formatted.txt",
    "SAS": "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/ldlc/ukb_sas_ldlc_common_snps_formatted.txt",
    "UGR": "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/ldlc/ugr_ldlc_common_snps_formatted.txt"
}

VCOR_PATHS = {
    "AFR": "/user/work/ac14629/MRC_network_project/data/UKB/LD_ref_dataset/AFR.multisusie.vcor",
    "AMR": "/user/work/ac14629/MRC_network_project/data/UKB/LD_ref_dataset/AMR.multisusie.vcor",
    "EAS": "/user/work/ac14629/MRC_network_project/data/UKB/LD_ref_dataset/EAS.multisusie.vcor",
    "EUR": "/user/work/ac14629/MRC_network_project/data/UKB/LD_ref_dataset/EUR.multisusie.vcor",
    "SAS": "/user/work/ac14629/MRC_network_project/data/UKB/LD_ref_dataset/SAS.multisusie.vcor",
    "UGR": "/user/work/ac14629/MRC_network_project/data/UKB/LD_ref_dataset/AFR.multisusie.vcor"
}

COL_MAP = {
    "AFR": {"chr": "chr", "pos": "pos", "Est": "beta", "Est.SE": "se"},
    "AMR": {"CHR": "chr", "POS": "pos", "BETA": "beta", "SE": "se"},
    "EAS": {"CHR": "chr", "BP": "pos", "BETA": "beta", "SE": "se"},
    "EUR": {"CHR": "chr", "BP": "pos", "BETA": "beta", "SE": "se"},
    "SAS": {"chr": "chr", "pos": "pos", "Est": "beta", "Est.SE": "se"},
    "UGR": {"chr": "chr", "pos": "pos", "Est": "beta", "Est.SE": "se"},
}

CHUNKSIZE = 1_000_000

# ----------------------------
# ARGUMENTS
# ----------------------------

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--snp", required=True)
    parser.add_argument("--window", type=int, default=1_000_000)
    parser.add_argument("--max_snps", type=int, default=5000)
    parser.add_argument("--outdir", default="results")

    return parser.parse_args()

# ----------------------------
# LOAD GWAS (chromosome only)
# ----------------------------

def load_gwas(pop, chrom):

    print(f"Loading {pop} GWAS chr{chrom}")

    usecols = list(COL_MAP[pop].keys()) + ["snp.id"]

    chr_col = [
        c for c in COL_MAP[pop]
        if COL_MAP[pop][c] == "chr"
    ][0]

    chunks = []

    for chunk in pd.read_csv(
        FILE_PATHS[pop],
        sep=r"\s+",
        usecols=usecols,
        chunksize=CHUNKSIZE
    ):

        chunk_chr = pd.to_numeric(
            chunk[chr_col],
            errors="coerce"
        )

        chunk = chunk[chunk_chr == chrom]

        if not chunk.empty:
            chunks.append(chunk)

    if len(chunks) == 0:
        raise ValueError(f"No chr{chrom} data found for {pop}")

    df = pd.concat(chunks, ignore_index=True)

    # harmonised SNP ID
    df["snp"] = df["snp.id"].astype(str)

    # rename columns
    df = df.rename(columns=COL_MAP[pop])

    df = df[["snp", "chr", "pos", "beta", "se"]]

    # check columns
    df["chr"] = pd.to_numeric(df["chr"], errors="coerce")
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    df["beta"] = pd.to_numeric(df["beta"], errors="coerce")
    df["se"] = pd.to_numeric(df["se"], errors="coerce")

    df = df.dropna()

    df = df.drop_duplicates("snp")

    print(f"{pop}: loaded {len(df)} chr{chrom} SNPs")

    return df

# ----------------------------
# FIND REGION
# ----------------------------

def get_region(target_snp, window):

    chrom, pos, *_ = target_snp.split(":")

    chrom = int(chrom)
    pos = int(pos)

    start = max(0, pos - window)
    end = pos + window

    return chrom, start, end

# ----------------------------
# FILTER REGION
# ----------------------------

def filter_region(df, chrom, start, end):

    return df[
        (df["chr"] == chrom) &
        (df["pos"] >= start) &
        (df["pos"] <= end)
    ]

# ----------------------------
# BUILD EFFECT AND SE LISTS
# ----------------------------

def build_effect_lists(data_dict, snps):

    b_list = []
    se_list = []

    for pop in POPULATIONS:

        df = (
            data_dict[pop]
            .set_index("snp")
            .loc[snps]
        )

        b_list.append(df["beta"].values)
        se_list.append(df["se"].values)

    return b_list, se_list

# ----------------------------
# BUILD LD MATRIX
# ----------------------------

def build_ld_matrix(pop, chrom, start, end, snps):

    print(f"Loading {pop} LD chr{chrom}")

    snp_set = frozenset(snps)

    usecols = [
        "#CHROM_A",
        "POS_A",
        "ID_A",
        "POS_B",
        "ID_B",
        "PHASED_R2"
    ]

    chunks = []

    for chunk in pd.read_csv(
        VCOR_PATHS[pop],
        sep="\t",
        usecols=usecols,
        chunksize=CHUNKSIZE
    ):

        chrom_a = pd.to_numeric(
            chunk["#CHROM_A"],
            errors="coerce"
        )

        pos_a = pd.to_numeric(
            chunk["POS_A"],
            errors="coerce"
        )

        pos_b = pd.to_numeric(
            chunk["POS_B"],
            errors="coerce"
        )

        # single combined mask
        mask = (
            (chrom_a == chrom) &
            (pos_a >= start) &
            (pos_a <= end) &
            (pos_b >= start) &
            (pos_b <= end)
        )

        chunk = chunk[mask]

        # DEBUG
        if not chunk.empty:
            print(f"\n{pop} rows after region filter:")
            print(chunk[["ID_A", "ID_B"]].head())
            print(f"Target SNP example: {list(snps)[:5]}")

        # SNP filter
        chunk = chunk[
            chunk["ID_A"].isin(snp_set) &
            chunk["ID_B"].isin(snp_set)
        ]

        if not chunk.empty:
            chunks.append(chunk)

    if len(chunks) == 0:

        raise ValueError(
            f"No LD data found for "
            f"{pop} chr{chrom}"
        )

    df = pd.concat(
        chunks,
        ignore_index=True
    )

    print(f"{pop}: retained {len(df)} LD pairs")

    snp_index = {
        s: i for i, s in enumerate(snps)
    }

    P = len(snps)

    R = np.zeros(
        (P, P),
        dtype=np.float32
    )

    i = df["ID_A"].map(snp_index).values
    j = df["ID_B"].map(snp_index).values

    r = np.sqrt(
        df["PHASED_R2"].values
    )

    R[i, j] = r
    R[j, i] = r

    np.fill_diagonal(R, 1.0)

    return R

# ----------------------------
# MAIN
# ----------------------------

def main():

    args = parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print(f"\nRunning SNP: {args.snp}")

    # determine chromosome immediately
    chrom, start, end = get_region(
        args.snp,
        args.window
    )

    print(
        f"Region: chr{chrom}:{start:,}-{end:,}"
    )

    # ----------------------------
    # LOAD GWAS
    # ----------------------------

    data_dict = {
        pop: load_gwas(pop, chrom)
        for pop in POPULATIONS
    }

    # ----------------------------
    # FILTER REGION
    # ----------------------------

    for pop in POPULATIONS:

        data_dict[pop] = filter_region(
            data_dict[pop],
            chrom,
            start,
            end
        )

        print(
            f"{pop}: "
            f"{len(data_dict[pop])} SNPs in region"
        )

    # ----------------------------
    # FIND COMMON SNPs
    # ----------------------------

    common_snps = set(
        data_dict["EUR"]["snp"]
    )

    for pop in POPULATIONS:
        common_snps &= set(
            data_dict[pop]["snp"]
        )

    common_snps = sorted(common_snps)

    print(
        f"Common SNPs across populations: "
        f"{len(common_snps)}"
    )

    if len(common_snps) == 0:
        print("No common SNPs — skipping")
        return

    if len(common_snps) > args.max_snps:
        print(
            f"Too many SNPs "
            f"({len(common_snps)}) — skipping"
        )
        return

    # ----------------------------
    # ALIGN TO COMMON SNPs
    # ----------------------------

    for pop in POPULATIONS:

        data_dict[pop] = (
            data_dict[pop]
            .set_index("snp")
            .loc[common_snps]
            .reset_index()
        )

    print(f"Final SNP count: {len(common_snps)}")

    # ----------------------------
    # BUILD EFFECTS
    # ----------------------------

    b_list, se_list = build_effect_lists(
        data_dict,
        common_snps
    )

    # ----------------------------
    # BUILD LD MATRICES
    # ----------------------------

    R_list = [
    build_ld_matrix(
        pop,
        chrom,
        start,
        end,
        common_snps
    )
    for pop in POPULATIONS
]

    # ----------------------------
    # POPULATION PARAMETERS
    # ----------------------------

    population_sizes_dict = {
        "AFR": 6225,
        "AMR": 34167,
        "EAS": 17687,
        "EUR": 467006,
        "SAS": 2629,
        "UGR": 6163
    }

    population_sizes = []

    for pop in POPULATIONS:

        if pop not in population_sizes_dict:
            raise ValueError(
                f"{pop} missing "
                f"from population_sizes_dict"
            )

        population_sizes.append(
            int(population_sizes_dict[pop])
        )

    # ----------------------------
    # LOAD varY
    # ----------------------------

    varY_df = pd.read_csv(
        "/user/work/ac14629/"
        "MRC_network_project/results/"
        "MAIN_ANALYSIS/ldlc/"
        "varY_estimates.csv"
    )

    varY_df["pop"] = (
        varY_df["pop"]
        .astype(str)
        .str.replace('"', '')
        .str.strip()
    )

    varY_df = varY_df.set_index("pop")

    varY_list = []
    missing = []

    for pop in POPULATIONS:

        if pop not in varY_df.index:
            missing.append(pop)
            continue

        varY_list.append(
            float(varY_df.loc[pop, "varY"])
        )

    if missing:
        raise ValueError(
            f"Missing populations in "
            f"varY_estimates.csv: {missing}"
        )

    print("\nPopulation parameters:")

    for pop, n, vy in zip(
        POPULATIONS,
        population_sizes,
        varY_list
    ):

        print(
            f"{pop}: "
            f"N={n}, "
            f"varY={vy:.4g}"
        )

    # ----------------------------
    # RUN MULTISUSIE
    # ----------------------------

    fit = MultiSuSiE.multisusie_rss(
        b_list=b_list,
        s_list=se_list,
        R_list=R_list,
        varY_list=varY_list,
        rho=0.75,
        population_sizes=population_sizes,
        L=10,
        low_memory_mode=True
    )

    # ----------------------------
    # SAVE OUTPUT
    # ----------------------------

    out = pd.DataFrame({
        "snp": common_snps,
        "pip": fit.pip
    })

    out_file = os.path.join(
        args.outdir,
        f"{args.snp}.pip.tsv"
    )

    out.to_csv(
        out_file,
        sep="\t",
        index=False
    )

    print(f"\nSaved: {out_file}")

# ----------------------------

if __name__ == "__main__":
    main()
