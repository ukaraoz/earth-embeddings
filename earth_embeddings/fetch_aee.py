import ee
import sys
import csv
import geemap
from datetime import datetime, timezone
from collections import defaultdict
from typing import List, Dict, Optional

EE_PROJECT = "embeddings-486622"

# Configuration AlphaEarth
COLLECTION_ID = 'GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL'
DEFAULT_YEAR = 2024
BAND_PREFIX = 'A'
DIMS = 64
BAND_LIST = [f'{BAND_PREFIX}{i:02d}' for i in range(DIMS)]
SCALE_M = 10.0


def log(msg):
    """Print message with timestamp."""
    ts = datetime.now().strftime('%H:%M:%S')
    print(f"[{ts}] {msg}", flush=True)


def init_gee():
    """Initialise Google Earth Engine."""
    log("Initialisation of GEE...")
    try:
        ee.Initialize()
        log("GEE initialized")
    except Exception:
        log("Authentification required")
        ee.Authenticate()
        ee.Initialize(project=EE_PROJECT)
        log("Authentification successful")


def load_locations(csv_path: str):
    """Load the latitude and longitudes for the locations"""
    log(f"Loading {csv_path}...")

    try:
        with open(csv_path, 'r', encoding='utf-8-sig', newline='') as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            fields = reader.fieldnames or []
    except FileNotFoundError:
        print(f"❌ File not found: {csv_path}")
        sys.exit(1)

    locations = []
    skipped_rows = []

    for r in rows:
        try:
            id_field = fields[0]
            id = r[id_field]
            lat = float(r["latitude"])
            lon = float(r["longitude"])

            locations.append({
                'id': id,
                'lat': lat,
                'lon': lon
            })
        except Exception as e:
            skipped_rows.append({**r, '_error': str(e)})
            continue

    if skipped_rows:
        skipped_path = csv_path.rsplit('.', 1)[0] + '_missinglatlong.csv'
        skipped_fields = list(skipped_rows[0].keys())
        with open(skipped_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=skipped_fields)
            writer.writeheader()
            writer.writerows(skipped_rows)
        log(f"⚠️  {len(skipped_rows)} skipped rows written to {skipped_path}")

    log(f"\n✅ {len(locations)} locations loaded ({len(skipped_rows)} ignored)")
    return locations


def get_year_mosaic(year: int):
    """Get AlphaEarth mosaic for the given year."""
    col = ee.ImageCollection(COLLECTION_ID)
    start = ee.Date.fromYMD(year, 1, 1)
    end = start.advance(1, 'year')

    filtered = col.filterDate(start, end)
    img = filtered.mosaic()

    # Fallback if year isn't available
    img = ee.Image(ee.Algorithms.If(
        filtered.size().gt(0),
        img,
        col.filterDate(
            start.advance(-3, 'year'),
            end.advance(3, 'year')
        ).mosaic()
    ))

    return ee.Image(img)


def sample_points_by_year(all_locations: List[dict], sample_year: int = DEFAULT_YEAR):
    """
    Add AlphaEarth to the locations.
    """
    results = {}

    # get the image
    img = get_year_mosaic(sample_year)

    # verify the bands
    band_names = img.bandNames().getInfo()
    log(f"  Available bands: {band_names[:5]}... ({len(band_names)} total)")

    # select the bands A00-A63
    img = img.select(BAND_LIST)
    # rename to f1-f64
    # img = img.rename([f'f{i}' for i in range(1, DIMS + 1)])

    # create FeatureCollection
    features = []
    for idx, pt in enumerate(all_locations):
        # print(f"idx={idx}, pt={pt}")
        features.append(ee.Feature(
            ee.Geometry.Point([pt['lon'], pt['lat']]),
            {'id': pt['id'], 'idx': idx}
        ))

    fc = ee.FeatureCollection(features)

    # can also use geemap
    # rgb_img = geemap.ee_to_numpy(img, region=fc)
    # print(rgb_img)

    # Adding embeddings
    log(f"   Ornamenting with embeddings...")
    sampled = img.sampleRegions(
        collection=fc,
        scale=SCALE_M,
        geometries=False,
        tileScale=1
    )
    # recover the results
    log(f"   Parsing the embeddings...")
    sampled_list = sampled.getInfo()

    if sampled_list and 'features' in sampled_list:
        n_results = len(sampled_list['features'])
        log(f"   ✅ {n_results} results received.")

        for feat in sampled_list['features']:
            props = feat.get('properties', {})
            idx = props.get('idx')
            id = props.get('id')
            if idx is not None:
                results[id] = props
            else:
                log(f"   ⚠️  No results")
    return results


def batch_sample_points_by_year(
    all_locations: List[dict],
    sample_year: int = DEFAULT_YEAR,
    batch_size: int = 500,
) -> Dict:
    """
    Split all_locations into batches and process each with sample_points_by_year,
    returning a merged dict keyed by biosample.
    """
    results = {}
    total = len(all_locations)
    n_batches = (total + batch_size - 1) // batch_size

    for batch_idx in range(n_batches):
        start = batch_idx * batch_size
        end = min(start + batch_size, total)
        batch = all_locations[start:end]

        log(f"Batch {batch_idx + 1}/{n_batches}: locations {start + 1}–{end} ({len(batch)} points)")
        try:
            batch_results = sample_points_by_year(batch, sample_year)
            results.update(batch_results)
            log(f"  Batch {batch_idx + 1} done — {len(batch_results)} results")
        except Exception as e:
            log(f"  ⚠️  Batch {batch_idx + 1} failed: {e}")

    log(f"\n✅ Batch processing complete: {len(results)}/{total} results returned.")
    return results


def write_output_csv(all_locations: List[dict], results: Dict, output_path: str):
    """
    Write output CSV

    Format: lat, lon, A00...A63
    """
    log(f"\Writing to {output_path}...")

    # Header
    header = ['id', 'lat', 'lon'] + BAND_LIST

    valid_count = 0
    invalid_count = 0
    skipped_rows = []
    skipped_path = output_path.rsplit('.', 1)[0] + '_skipped.csv'

    with open(output_path, 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)

        for idx, pt in enumerate(all_locations):
            props = results.get(pt['id'], {})

            # extract the features
            features = []
            row_errors = []

            for i in range(DIMS):
                band = f'{BAND_PREFIX}{i:02d}'
                val = props.get(band)
                try:
                    features.append(float(val))
                except Exception as e:
                    row_errors.append(f'{band}: {e}')

            # verify if we have valid results
            # f==f is F for NaN
            n_valid = sum(1 for f in features if f == f)

            if n_valid > 0:
                valid_count += 1
            else:
                invalid_count += 1

            # write the lign
            row = [
                pt['id'],
                f"{pt['lat']:.6f}",
                f"{pt['lon']:.6f}"
            ] + [f"{f:.15f}" if f == f else "" for f in features]

            writer.writerow(row)

            if row_errors:
                skipped_rows.append({
                    'id': pt['id'],
                    'lat': pt['lat'],
                    'lon': pt['lon'],
                    '_errors': '; '.join(row_errors),
                })

    if skipped_rows:
        with open(skipped_path, 'w', encoding='utf-8', newline='') as sf:
            skipped_writer = csv.DictWriter(
                sf, fieldnames=['id', 'lat', 'lon', '_errors'])
            skipped_writer.writeheader()
            skipped_writer.writerows(skipped_rows)
        log(f"⚠️  {len(skipped_rows)} rows with missing/invalid bands written to {skipped_path}")

    log(f"\n✅ Finished.")
    log(f"   Valid points: {valid_count}/{len(all_locations)}")
    log(f"   Missing data points: {invalid_count}/{len(all_locations)}")


def main():
    """Main function"""
    print("\n" + "="*80)
    print("  Generating CSV for locations.csv - AlphaEarth embedding vectors")
    print("="*80)

    # Configuration
    # LOCATIONS_CSV = 'data/locations-small.csv'
    # OUTPUT_CSV = 'data/locations_with_vectors.csv'
    LOCATIONS_CSV = 'data/USGS_soil_geochemfeatures_cols1_2_18_19.csv'
    OUTPUT_CSV = 'data/USGS_soil_geochemfeatures_cols1_2_18_19_withvectors.csv'

    # 1. Initialisation
    init_gee()

    # 2. Charger les événements
    locations = load_locations(LOCATIONS_CSV)

    if not locations:
        print("❌ Missing valid locations.")
        sys.exit(1)

    log(f"\n📊 Total: {len(locations)} locations")

    # results = sample_points_by_year(locations)
    results = batch_sample_points_by_year(locations, batch_size=2000)
    # 4. Écrire le CSV
    write_output_csv(locations, results, OUTPUT_CSV)


if __name__ == "__main__":
    main()
