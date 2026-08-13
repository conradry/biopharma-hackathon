#!/usr/bin/env python3
"""
toxin_pd_correlation.py

Build a county x year panel linking environmental toxin exposure to
Parkinson's disease (PD) incidence, and compute lagged correlations
between them -- i.e. hypothesis generation for
"does exposure to chemical X predict a later rise in PD incidence?"

DATA SOURCES
------------
1. Pesticide exposure (automated):
   USGS Pesticide National Synthesis Project, county-level annual use
   estimates, 1992-2012 (tab-delimited files, no API key needed):
   https://water.usgs.gov/nawqa/pnsp/usage/maps/county-level/PesticideUseEstimates/EPest.county.estimates.{year}.txt

2. Air pollution (automated):
   CDC Environmental Public Health Tracking Network, county-level PM2.5,
   via the Socrata Open Data API (data.cdc.gov), dataset dqwm-pbi7
   ("Daily County-Level PM2.5 Concentrations, 2001-2019"):
   https://data.cdc.gov/resource/dqwm-pbi7.json

3. Parkinson's disease incidence (MANUAL -- no public API):
   IHME Global Burden of Disease (GBD) results require a free account
   and export from the GBD Results Tool:
   https://vizhub.healthdata.org/gbd-results/
   Query: Cause = "Parkinson's disease", Measure = "Incidence",
          Metric = "Rate", Location = US states/counties if available,
          Year = as many years as possible, Sex = Both, Age = Age-standardized.
   Export as CSV and pass its path via --pd-csv. Required columns
   (rename yours to match, or pass --pd-* column-name overrides):
       county_fips, year, pd_incidence_rate

USAGE
-----
# Step 1: pull pesticide + PM2.5 data for a set of years
python3 toxin_pd_correlation.py fetch --years 2000 2012 --out panel_exposure.csv

# Step 2: merge with your manually-exported GBD PD incidence CSV and
# compute lagged correlations
python3 toxin_pd_correlation.py correlate \
    --exposure panel_exposure.csv \
    --pd-csv gbd_pd_incidence.csv \
    --lags 0 5 10 15 20 \
    --out toxin_pd_correlations.csv
"""

import argparse
import sys
import time
import numpy as np
import pandas as pd
import requests
from scipy.stats import pearsonr

USGS_PESTICIDE_URL = (
    "https://water.usgs.gov/nawqa/pnsp/usage/maps/county-level/"
    "PesticideUseEstimates/EPest.county.estimates.{year}.txt"
)
CDC_PM25_API = "https://data.cdc.gov/resource/dqwm-pbi7.json"


# ---------------------------------------------------------------------
# Step 1a: USGS pesticide data
# ---------------------------------------------------------------------

def fetch_pesticide_year(year, compounds=None):
    """Download and parse one year of USGS county-level pesticide-use estimates."""
    url = USGS_PESTICIDE_URL.format(year=year)
    print(f"  fetching pesticide data for {year} ...")
    resp = requests.get(url, timeout=60)
    if resp.status_code != 200:
        print(f"  ! no USGS file for {year} (status {resp.status_code}); "
              f"years 1992-2012 are available via this URL pattern, "
              f"2013+ requires the ScienceBase downloads instead", file=sys.stderr)
        return pd.DataFrame()

    from io import StringIO
    df = pd.read_csv(StringIO(resp.text), sep="\t", dtype=str)
    df.columns = [c.strip() for c in df.columns]
    df["FIPS"] = df["STATE_FIPS_CODE"].str.zfill(2) + df["COUNTY_FIPS_CODE"].str.zfill(3)
    df["YEAR"] = df["YEAR"].astype(int)
    # average of low/high estimate as a single exposure figure, in kg
    df["EPEST_LOW_KG"] = pd.to_numeric(df["EPEST_LOW_KG"], errors="coerce")
    df["EPEST_HIGH_KG"] = pd.to_numeric(df["EPEST_HIGH_KG"], errors="coerce")
    df["EPEST_KG"] = df[["EPEST_LOW_KG", "EPEST_HIGH_KG"]].mean(axis=1)

    if compounds:
        df = df[df["COMPOUND"].str.lower().isin([c.lower() for c in compounds])]

    return df[["FIPS", "YEAR", "COMPOUND", "EPEST_KG"]]


def fetch_pesticides(years, compounds=None, pause=0.5):
    frames = []
    for year in years:
        frame = fetch_pesticide_year(year, compounds=compounds)
        if not frame.empty:
            frames.append(frame)
        time.sleep(pause)
    if not frames:
        return pd.DataFrame(columns=["FIPS", "YEAR", "COMPOUND", "EPEST_KG"])
    combined = pd.concat(frames, ignore_index=True)
    # pivot: one column per compound
    wide = combined.pivot_table(
        index=["FIPS", "YEAR"], columns="COMPOUND", values="EPEST_KG", aggfunc="sum"
    ).reset_index()
    wide.columns = [
        c if c in ("FIPS", "YEAR") else f"pesticide_{c.replace(' ', '_')}"
        for c in wide.columns
    ]
    return wide


# ---------------------------------------------------------------------
# Step 1b: CDC PM2.5 data (Socrata API)
# ---------------------------------------------------------------------

def fetch_pm25(years, pause=0.3):
    """
    Pull county-level daily PM2.5 from CDC's Socrata API and aggregate to
    annual county means. Column names are auto-detected defensively since
    Socrata field names can change between dataset revisions -- inspect
    the first response if this needs adjusting.
    """
    print("  fetching PM2.5 data from CDC (data.cdc.gov) ...")
    frames = []
    for year in years:
        params = {
            "$where": f"date_extract_y(date)={year}",
            "$select": "countyfips,date,pm25_max_pred",
            "$limit": 50000,
        }
        try:
            resp = requests.get(CDC_PM25_API, params=params, timeout=60)
            resp.raise_for_status()
            data = resp.json()
        except requests.RequestException as err:
            print(f"  ! PM2.5 fetch failed for {year}: {err}", file=sys.stderr)
            continue
        if not data:
            continue
        df = pd.DataFrame(data)
        # defensive column detection
        fips_col = next((c for c in df.columns if "fips" in c.lower()), None)
        pm_col = next((c for c in df.columns if "pm25" in c.lower()), None)
        if not fips_col or not pm_col:
            print(f"  ! unexpected PM2.5 schema for {year}: {list(df.columns)}", file=sys.stderr)
            continue
        df[pm_col] = pd.to_numeric(df[pm_col], errors="coerce")
        annual = df.groupby(fips_col)[pm_col].mean().reset_index()
        annual.columns = ["FIPS", "pm25_annual_mean"]
        annual["FIPS"] = annual["FIPS"].astype(str).str.zfill(5)
        annual["YEAR"] = year
        frames.append(annual)
        time.sleep(pause)
    if not frames:
        return pd.DataFrame(columns=["FIPS", "YEAR", "pm25_annual_mean"])
    return pd.concat(frames, ignore_index=True)


# ---------------------------------------------------------------------
# Step 2: merge with PD incidence and compute lagged correlations
# ---------------------------------------------------------------------

def load_pd_incidence(path, fips_col="county_fips", year_col="year", rate_col="pd_incidence_rate"):
    df = pd.read_csv(path, dtype={fips_col: str})
    df = df.rename(columns={fips_col: "FIPS", year_col: "YEAR", rate_col: "PD_INCIDENCE"})
    df["FIPS"] = df["FIPS"].str.zfill(5)
    df["YEAR"] = df["YEAR"].astype(int)
    return df[["FIPS", "YEAR", "PD_INCIDENCE"]]


def lagged_correlation(panel, exposure_col, lags):
    """
    For a given exposure column, compute Pearson correlation between
    exposure at year t and PD incidence at year t+lag, pooling all
    counties. Returns a DataFrame of lag, r, p, n.
    """
    rows = []
    for lag in lags:
        shifted = panel[["FIPS", "YEAR", exposure_col]].copy()
        shifted["YEAR"] = shifted["YEAR"] + lag  # align exposure year t to PD year t+lag
        merged = shifted.merge(
            panel[["FIPS", "YEAR", "PD_INCIDENCE"]], on=["FIPS", "YEAR"], how="inner"
        ).dropna(subset=[exposure_col, "PD_INCIDENCE"])
        if len(merged) < 10:
            rows.append({"lag": lag, "r": np.nan, "p": np.nan, "n": len(merged)})
            continue
        r, p = pearsonr(merged[exposure_col], merged["PD_INCIDENCE"])
        rows.append({"lag": lag, "r": r, "p": p, "n": len(merged)})
    return pd.DataFrame(rows)


def run_correlations(panel, lags):
    exposure_cols = [
        c for c in panel.columns
        if c.startswith("pesticide_") or c == "pm25_annual_mean"
    ]
    results = []
    for col in exposure_cols:
        lag_df = lagged_correlation(panel, col, lags)
        lag_df["exposure"] = col
        results.append(lag_df)
    all_results = pd.concat(results, ignore_index=True)
    return all_results.sort_values("r", key=lambda s: s.abs(), ascending=False)


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    p_fetch = sub.add_parser("fetch", help="pull pesticide + PM2.5 data into a county x year panel")
    p_fetch.add_argument("--years", nargs=2, type=int, metavar=("START", "END"), required=True)
    p_fetch.add_argument("--compounds", nargs="+", default=None,
                          help="restrict to specific pesticide compounds, e.g. paraquat rotenone")
    p_fetch.add_argument("--out", default="panel_exposure.csv")

    p_corr = sub.add_parser("correlate", help="merge exposure panel with PD incidence and compute lagged correlations")
    p_corr.add_argument("--exposure", required=True, help="CSV from the 'fetch' step")
    p_corr.add_argument("--pd-csv", required=True, help="manually exported GBD PD incidence CSV")
    p_corr.add_argument("--pd-fips-col", default="county_fips")
    p_corr.add_argument("--pd-year-col", default="year")
    p_corr.add_argument("--pd-rate-col", default="pd_incidence_rate")
    p_corr.add_argument("--lags", nargs="+", type=int, default=[0, 5, 10, 15, 20])
    p_corr.add_argument("--out", default="toxin_pd_correlations.csv")

    args = parser.parse_args()

    if args.command == "fetch":
        years = list(range(args.years[0], args.years[1] + 1))
        pesticide_df = fetch_pesticides(years, compounds=args.compounds)
        pm25_df = fetch_pm25(years)
        panel = pesticide_df.merge(pm25_df, on=["FIPS", "YEAR"], how="outer")
        panel.to_csv(args.out, index=False)
        print(f"\nSaved exposure panel ({len(panel)} rows, {len(panel.columns)} columns) to {args.out}")
        print("NOTE: this covers pesticide use (1992-2012 via USGS static files) and PM2.5.")
        print("Next: export PD incidence from https://vizhub.healthdata.org/gbd-results/ "
              "and run the 'correlate' command.")

    elif args.command == "correlate":
        panel = pd.read_csv(args.exposure, dtype={"FIPS": str})
        panel["FIPS"] = panel["FIPS"].str.zfill(5)
        pd_df = load_pd_incidence(args.pd_csv, args.pd_fips_col, args.pd_year_col, args.pd_rate_col)
        panel = panel.merge(pd_df, on=["FIPS", "YEAR"], how="left")

        if panel["PD_INCIDENCE"].notna().sum() == 0:
            print("! No matching FIPS/YEAR rows between exposure panel and PD incidence file. "
                  "Check that FIPS codes and years line up.", file=sys.stderr)
            sys.exit(1)

        results = run_correlations(panel, args.lags)
        results.to_csv(args.out, index=False)
        print(f"\nSaved {len(results)} lag x exposure correlation rows to {args.out}")
        print("\nTop 15 by |r| (hypothesis-generating only -- NOT causal evidence):")
        print(results.head(15).to_string(index=False))
        print(
            "\nReminder: county-level correlations are confounded by age structure, "
            "rurality, smoking, healthcare access, and other co-occurring exposures. "
            "Treat strong hits as candidates for follow-up, not conclusions."
        )


if __name__ == "__main__":
    main()
