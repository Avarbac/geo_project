# Геометрические алгоритмы на данных OpenStreetMap

## Данные
- `data/points.json` — 818 остановок `highway=bus_stop` из Overpass API, отметка данных `2026-10-02T20:21:34Z`.
- `data/district.geojson` — Кировский район Уфы, OSM relation 5523570 (GeoJSON извлечён из приложенного `code(3).txt`).

## Установка
Python 3.10+.
python -m venv .venv
pip install -r requirements.txt
```

## Запуск из корня проекта
```bash
python geo.py project data/points.json data/district.geojson
python geo.py hull results/projected/points.txt
python geo.py closest results/projected/points.txt
python geo.py closest results/projected/points.txt --brute --output results/closest_brute.json
python geo.py inside data/points.json data/district.geojson
python geo.py plot data/points.json data/district.geojson
python validate.py
python -m unittest discover -s tests -v
```

## Алгоритмы
- Проекция: локальная равнопромежуточная, с началом в центре охватывающего прямоугольника данных, округление до целого метра.
- Выпуклая оболочка: монотонная цепочка, `O(n log n)`.
- Площадь: формула шнурования.
- Принадлежность району: пересечение луча с кольцом; для Polygon учитываются внутренние кольца-отверстия, для MultiPolygon — все компоненты.
- Ближайшая пара: divide-and-conquer `O(n log n)`; `--brute` — эталонный перебор `O(n²)`.
- Основная программа не зависит от Shapely; `validate.py` применяет Shapely/GEOS только для независимого сравнения.

## Форматы промежуточных файлов
`results/projected/points.txt`: `OSM_ID x_m y_m`, одна точка в строке.
`results/projected/district.txt`: координаты колец по строкам; пустая строка разделяет кольца и компоненты. Внешние кольца нормализуются против часовой стрелки, отверстия — по часовой.
