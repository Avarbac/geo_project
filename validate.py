"""Независимая проверка результатов через Shapely (не используется основными алгоритмами)."""
import json,sys
from pathlib import Path
from shapely.geometry import Point, shape
from geo import read_points, read_geojson, projection_origin, project_lonlat, convex_hull, closest_pair, load_projected

def main():
    points_path=sys.argv[1] if len(sys.argv)>1 else 'data/points.json' 
    district_path=sys.argv[2] if len(sys.argv)>2 else 'data/district.geojson'
    pts=read_points(points_path); polygons=read_geojson(district_path)
    phi0,lam0=projection_origin(pts,polygons)
    coords=[project_lonlat(p['lon'],p['lat'],phi0,lam0) for p in pts]
    geo=json.loads(Path(district_path).read_text(encoding='utf8')); features=geo.get('features',[geo])
    # Build projected polygon from same projection and compare against independent GEOS predicates.
    from shapely.geometry import Polygon, MultiPolygon
    sh_polys=[]
    for poly in polygons:
        rings=[]
        for ring in poly:
            rr=[project_lonlat(lon,lat,phi0,lam0) for lon,lat,*_ in ring]
            rings.append(rr)
        sh_polys.append(Polygon(rings[0],rings[1:]))
    shp=MultiPolygon(sh_polys)
    from geo import point_in_multipolygon
    projected=[]
    for poly in polygons:
        projected.append([[project_lonlat(lon,lat,phi0,lam0) for lon,lat,*_ in ring] for ring in poly])
    manual=[point_in_multipolygon(projected,p) for p in coords]
    independent=[shp.covers(Point(p)) for p in coords]
    mismatches=[i for i,(a,b) in enumerate(zip(manual,independent)) if a!=b]
    hull=convex_hull(coords);sh_hull=__import__('shapely.geometry',fromlist=['MultiPoint']).MultiPoint(coords).convex_hull
    pair,d,checks=closest_pair(coords)
    print(json.dumps({'points':len(pts),'inside_manual':sum(manual),'inside_shapely':sum(independent),'inside_mismatches':len(mismatches),'mismatch_indices':mismatches[:20],'manual_hull_vertices':len(hull),'shapely_hull_vertices':len(list(sh_hull.exterior.coords))-1,'manual_hull_area_m2':abs(__import__('geo').signed_area2(hull))/2,'shapely_hull_area_m2':sh_hull.area,'closest_pair_distance_m':d,'closest_pair_checks':checks},ensure_ascii=False,indent=2))
if __name__=='__main__':main()
