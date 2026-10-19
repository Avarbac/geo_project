"""Учебная CLI-утилита геометрии для точек OSM. Без Shapely в основном алгоритме."""
from __future__ import annotations
import argparse, json, math, time, csv
from pathlib import Path
from typing import Iterable

R_EARTH = 6_371_000.0

def read_points(path):
    data=json.loads(Path(path).read_text(encoding='utf-8'))
    if isinstance(data,dict) and 'elements' in data:
        rows=[]
        for e in data['elements']:
            if 'lat' in e and 'lon' in e:
                tags=e.get('tags',{})
                rows.append({'id':str(e.get('id','')), 'lat':float(e['lat']), 'lon':float(e['lon']), 'name':tags.get('name',''), 'tags':tags})
            elif 'center' in e:
                rows.append({'id':str(e.get('id','')), 'lat':float(e['center']['lat']), 'lon':float(e['center']['lon']), 'name':e.get('tags',{}).get('name',''), 'tags':e.get('tags',{})})
        return rows
    if isinstance(data,list):
        return [{'id':str(p.get('id',i)), 'lat':float(p['lat']), 'lon':float(p['lon']), 'name':p.get('name',''), 'tags':p.get('tags',{})} for i,p in enumerate(data)]
    raise ValueError('Ожидается JSON Overpass с elements или список точек {id,lat,lon}')

def read_geojson(path):
    obj=json.loads(Path(path).read_text(encoding='utf-8'))
    if obj.get('type')=='FeatureCollection':
        features=obj.get('features',[])
        if not features: raise ValueError('В GeoJSON нет features')
        geoms=[f['geometry'] for f in features if f.get('geometry')]
    elif obj.get('type')=='Feature': geoms=[obj['geometry']]
    else: geoms=[obj]
    polygons=[]
    for g in geoms:
        if g['type']=='Polygon': polygons.append(g['coordinates'])
        elif g['type']=='MultiPolygon': polygons.extend(g['coordinates'])
        else: raise ValueError('Поддерживаются только Polygon и MultiPolygon')
    return polygons

def projection_origin(points, polygons):
    lats=[p['lat'] for p in points]; lons=[p['lon'] for p in points]
    for poly in polygons:
        for ring in poly:
            for lon,lat,*_ in ring: lats.append(lat); lons.append(lon)
    return (min(lats)+max(lats))/2, (min(lons)+max(lons))/2

def project_lonlat(lon,lat,phi0,lam0):
    return (round(R_EARTH*math.radians(lon-lam0)*math.cos(math.radians(phi0))), round(R_EARTH*math.radians(lat-phi0)))

def orient(a,b,c): return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
def dist2(a,b): return (a[0]-b[0])**2+(a[1]-b[1])**2

def convex_hull(points):
    pts=sorted(set(tuple(p) for p in points))
    if len(pts)<=1:return pts
    lo=[]
    for p in pts:
        while len(lo)>=2 and orient(lo[-2],lo[-1],p)<=0:lo.pop()
        lo.append(p)
    hi=[]
    for p in reversed(pts):
        while len(hi)>=2 and orient(hi[-2],hi[-1],p)<=0:hi.pop()
        hi.append(p)
    return lo[:-1]+hi[:-1]

def signed_area2(ring):
    return sum(ring[i][0]*ring[(i+1)%len(ring)][1]-ring[(i+1)%len(ring)][0]*ring[i][1] for i in range(len(ring))) if len(ring)>=3 else 0

def point_in_ring(p,ring):
    x,y=p; inside=False
    if len(ring)<3:return False
    for i in range(len(ring)):
        x1,y1=ring[i]; x2,y2=ring[(i+1)%len(ring)]
        if orient((x1,y1),(x2,y2),(x,y))==0 and min(x1,x2)<=x<=max(x1,x2) and min(y1,y2)<=y<=max(y1,y2): return True
        if (y1>y)!=(y2>y):
            xcross=x1+(y-y1)*(x2-x1)/(y2-y1)
            if xcross>=x: inside=not inside
    return inside

def point_in_polygon(polygon,p):
    if not polygon or not point_in_ring(p,polygon[0]):return False
    return not any(point_in_ring(p,hole) for hole in polygon[1:])
def point_in_multipolygon(polygons,p): return any(point_in_polygon(poly,p) for poly in polygons)

def closest_slow(points):
    n=len(points); best=None; pair=None; checks=0
    for i in range(n):
        for j in range(i+1,n):
            checks+=1; d=dist2(points[i],points[j])
            if best is None or d<best:best=d;pair=(i,j)
    return pair, (math.sqrt(best) if best is not None else math.inf), checks

def closest_pair(points):
    """Divide and conquer O(n log n); returns original point indices, distance, distance checks."""
    if len(points)<2:return None,math.inf,0
    px=sorted([(p[0],p[1],i) for i,p in enumerate(points)])
    py=sorted(px,key=lambda z:(z[1],z[0],z[2])); checks=0
    def d2(a,b):
        nonlocal checks; checks+=1
        return (a[0]-b[0])**2+(a[1]-b[1])**2
    def rec(xs,ys):
        n=len(xs)
        if n<=3:
            best=float('inf'); pair=None
            for i in range(n):
                for j in range(i+1,n):
                    dd=d2(xs[i],xs[j])
                    if dd<best:best=dd;pair=(xs[i][2],xs[j][2])
            return best,pair
        mid=n//2; left=xs[:mid]; right=xs[mid:]; splitx=right[0][0]
        li={z[2] for z in left}; ly=[z for z in ys if z[2] in li]; ry=[z for z in ys if z[2] not in li]
        dl,pl=rec(left,ly); dr,pr=rec(right,ry)
        best,pair=(dl,pl) if dl<=dr else (dr,pr)
        strip=[z for z in ys if (z[0]-splitx)**2<best]
        for i in range(len(strip)):
            j=i+1
            while j<len(strip) and (strip[j][1]-strip[i][1])**2<best:
                dd=d2(strip[i],strip[j])
                if dd<best:best=dd;pair=(strip[i][2],strip[j][2])
                j+=1
        return best,pair
    best,pair=rec(px,py)
    return pair,math.sqrt(best),checks

def save_projected(points, polygons, outdir):
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    phi0,lam0=projection_origin(points,polygons)
    coords=[project_lonlat(p['lon'],p['lat'],phi0,lam0) for p in points]
    with (out/'points.txt').open('w',encoding='utf-8') as f:
        for p,(x,y) in zip(points,coords):f.write(f"{p['id']} {x} {y}\n")
    projected=[]
    for poly in polygons:
        rings=[]
        for ri,ring in enumerate(poly):
            rr=[project_lonlat(lon,lat,phi0,lam0) for lon,lat,*_ in ring]
            if rr and rr[0]==rr[-1]:rr=rr[:-1]
            should_positive=(ri==0)
            if (signed_area2(rr)>0)!=should_positive:rr=list(reversed(rr))
            rings.append(rr)
        projected.append(rings)
    with (out/'district.txt').open('w',encoding='utf-8') as f:
        for pi,poly in enumerate(projected):
            if pi:f.write('\n')
            for ri,ring in enumerate(poly):
                if ri:f.write('\n')
                for x,y in ring:f.write(f'{x} {y}\n')
    meta={'phi0':phi0,'lambda0':lam0,'earth_radius_m':R_EARTH,'points_count':len(points),'polygons_count':len(polygons),'projection':'local equirectangular; x=R*(lon-lambda0)*cos(phi0), y=R*(lat-phi0); rounded to integer metres'}
    (out/'projection.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding='utf-8')
    return coords,projected,meta

def write_json(path,obj):Path(path).write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf-8')
def area_m2(ring):return abs(signed_area2(ring))/2

def command_project(args):
    pts=read_points(args.points); polys=read_geojson(args.district); coords,rings,meta=save_projected(pts,polys,args.outdir)
    print(f"Проекция завершена: точек {len(pts)}, полигонов {len(polys)}, начало координат φ₀={meta['phi0']:.8f}, λ₀={meta['lambda0']:.8f}")
    print(f"Файлы: {args.outdir}/points.txt, district.txt, projection.json")

def load_projected(path):
    ids=[];points=[]
    for line in Path(path).read_text(encoding='utf-8').splitlines():
        s=line.split()
        if len(s)>=3:ids.append(s[0]);points.append((int(s[1]),int(s[2])))
    return ids,points

def command_hull(args):
    ids,pts=load_projected(args.points); hull=convex_hull(pts); area=area_m2(hull)/1e6
    Path(args.output).parent.mkdir(parents=True,exist_ok=True)
    Path(args.output).write_text(''.join(f'{x} {y}\n' for x,y in hull),encoding='utf-8')
    left=min(pts,key=lambda p:p[0]) if pts else None
    obj={'points':len(pts),'unique_coordinates':len(set(pts)),'hull_vertices':len(hull),'area_km2':area,'westmost_point':left,'hull':hull}
    write_json(args.output+'.json',obj)
    print(f"Точек: {len(pts)}; уникальных координат: {len(set(pts))}; вершин оболочки: {len(hull)}; площадь: {area:.6f} км²")
    if left:print(f"Самая западная точка входа: {left}; входит в оболочку: {left in hull}")
    print(f"Оболочка записана в {args.output}")

def command_closest(args):
    ids,pts=load_projected(args.points); fn=closest_slow if args.brute else closest_pair
    t=time.perf_counter(); pair,d,checks=fn(pts); elapsed=time.perf_counter()-t
    if pair is None: print('Нужно минимум две точки');return
    i,j=pair
    result={'algorithm':'brute force O(n²)' if args.brute else 'divide and conquer O(n log n)','points':len(pts),'index_pair':[i,j],'id_pair':[ids[i],ids[j]],'point_pair':[pts[i],pts[j]],'distance_m':d,'distance_checks':checks,'elapsed_seconds':elapsed}
    write_json(args.output,result)
    print(f"Алгоритм: {result['algorithm']}; ближайшая пара: OSM {ids[i]} {pts[i]} — OSM {ids[j]} {pts[j]}; расстояние {d:.3f} м; проверок расстояния {checks}; время {elapsed:.4f} с")
    print(f"Результат: {args.output}")

def command_inside(args):
    pts=read_points(args.points); polys=read_geojson(args.district); coords,projected,meta=save_projected(pts,polys,args.outdir)
    inside=[]
    for p,xy in zip(pts,coords):
        if point_in_multipolygon(projected,xy):inside.append({**{k:p[k] for k in ('id','name','lat','lon')},'x':xy[0],'y':xy[1]})
    write_json(args.output,{'input_points':len(pts),'inside_count':len(inside),'outside_count':len(pts)-len(inside),'inside_points':inside})
    print(f"Внутри границы: {len(inside)} из {len(pts)}; снаружи: {len(pts)-len(inside)}; результат: {args.output}")

def command_plot(args):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    pts_raw=read_points(args.points);polys=read_geojson(args.district); coords,projected,meta=save_projected(pts_raw,polys,args.outdir)
    hull=convex_hull(coords);pair,d,checks=closest_pair(coords)
    inside=[point_in_multipolygon(projected,p) for p in coords]
    fig,ax=plt.subplots(figsize=(10,9));
    ax.scatter([p[0] for p in coords],[p[1] for p in coords],s=9,c=['#2563eb' if v else '#9ca3af' for v in inside],alpha=.7,label='Остановки: внутри/снаружи')
    for poly in projected:
        for ri,ring in enumerate(poly):
            loop=ring+[ring[0]] if ring else []
            ax.plot([p[0] for p in loop],[p[1] for p in loop],color='#15803d',linewidth=1.2,label='Граница района' if ri==0 else None)
    if hull:
        loop=hull+[hull[0]];ax.plot([p[0] for p in loop],[p[1] for p in loop],color='#dc2626',linestyle='--',linewidth=1.5,label='Выпуклая оболочка всех точек')
    if pair:
        a,b=coords[pair[0]],coords[pair[1]];ax.plot([a[0],b[0]],[a[1],b[1]],color='#9333ea',linewidth=2.5,label=f'Ближайшая пара ({d:.1f} м)');ax.scatter([a[0],b[0]],[a[1],b[1]],c='#9333ea',s=35)
    ax.set_aspect('equal',adjustable='datalim');ax.set_xlabel('x, м');ax.set_ylabel('y, м');ax.set_title('Остановки OSM и Кировский район Уфы');ax.grid(alpha=.2);ax.legend(loc='best',fontsize=8);fig.tight_layout()
    out=Path(args.output);out.parent.mkdir(parents=True,exist_ok=True);fig.savefig(out,dpi=180);plt.close(fig)
    print(f"График сохранён: {out}")

def build_parser():
    p=argparse.ArgumentParser(description='Геометрические алгоритмы для точек OSM');sub=p.add_subparsers(dest='command',required=True)
    a=sub.add_parser('project',help='проецировать точки и район в целые метры');a.add_argument('points');a.add_argument('district');a.add_argument('--outdir',default='results/projected');a.set_defaults(func=command_project)
    a=sub.add_parser('hull',help='выпуклая оболочка');a.add_argument('points',help='файл points.txt после project');a.add_argument('--output',default='results/hull.txt');a.set_defaults(func=command_hull)
    a=sub.add_parser('closest',help='ближайшая пара');a.add_argument('points',help='файл points.txt после project');a.add_argument('--brute',action='store_true',help='полный перебор O(n²) вместо divide-and-conquer');a.add_argument('--output',default='results/closest.json');a.set_defaults(func=command_closest)
    a=sub.add_parser('inside',help='точки внутри района');a.add_argument('points');a.add_argument('district');a.add_argument('--outdir',default='results/projected');a.add_argument('--output',default='results/inside.json');a.set_defaults(func=command_inside)
    a=sub.add_parser('plot',help='визуализация точек, района, оболочки и ближайшей пары');a.add_argument('points');a.add_argument('district');a.add_argument('--outdir',default='results/projected');a.add_argument('--output',default='results/map.png');a.set_defaults(func=command_plot)
    return p
if __name__=='__main__':
    args=build_parser().parse_args();args.func(args)
