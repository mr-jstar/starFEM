#!/usr/bin/env python

import argparse
import math
from typing import List, Tuple, Optional

import gmsh

def parse_box(s: str):
    """
    Format:
        "[x_min:x_max]x[y_min:y_max]x[z_min:z_max]"

    Przykład:
        "[0:100]x[0:100]x[-20:0]"
    """
    s = s.strip()

    try:
        parts = s.split("x")
        if len(parts) != 3:
            raise ValueError

        ranges = []
        for part in parts:
            part = part.strip()
            if not (part.startswith("[") and part.endswith("]")):
                raise ValueError

            part = part[1:-1]
            vmin, vmax = part.split(":")
            ranges.append((float(vmin), float(vmax)))

        (x_min, x_max), (y_min, y_max), (z_min, z_max) = ranges

        if x_max <= x_min or y_max <= y_min or z_max <= z_min:
            raise ValueError

        return x_min, x_max, y_min, y_max, z_min, z_max

    except Exception:
        raise argparse.ArgumentTypeError(
            'Niepoprawny format --box. Użyj np. "[0:100]x[0:100]x[-20:0]"'
        )

def auto_positions(
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
    rods: List[Tuple[float, float]],
):
    n = len(rods)

    if n == 0:
        return []

    width = x_max - x_min
    height = y_max - y_min

    nx = math.ceil(math.sqrt(n * width / height))
    ny = math.ceil(n / nx)

    dx = width / (nx + 1)
    dy = height / (ny + 1)

    xy = []

    for j in range(ny):
        for i in range(nx):
            if len(xy) >= n:
                break

            x = x_min + (i + 1) * dx
            y = y_min + (j + 1) * dy

            xy.append((x, y))

    return xy

def parse_rods(s: str):
    """
    Format:
        "0.016:0.8,0.02:1.2,0.012:0.6"
    czyli:
        d:l,d:l,d:l
    """
    rods = []
    for item in s.split(","):
        d, l = item.split(":")
        rods.append((float(d), float(l)))
    return rods


def parse_xy(s: str):
    """
    Format:
        "20:20,50:50,80:70"
    czyli:
        x:y,x:y,x:y
    """
    xy = []
    for item in s.split(","):
        x, y = item.split(":")
        xy.append((float(x), float(y)))
    return xy


def add_physical_surface(name: str, tags: List[int], phys_id: int):
    tags = sorted(set(tags))
    if tags:
        gmsh.model.addPhysicalGroup(2, tags, phys_id)
        gmsh.model.setPhysicalName(2, phys_id, name)


def generate_mesh(
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
    z_min: float,
    z_max: float,
    rods: List[Tuple[float, float]],
    xy: Optional[List[Tuple[float, float]]],
    output: str,
    lc_ground: float,
    lc_rods: float,
    msh_version: float = 2.2,
):
    if xy is None:
        xy = auto_positions(x_min, x_max, y_min, y_max, rods)

    if len(xy) != len(rods):
        raise ValueError("Liczba pozycji xy musi być taka sama jak liczba prętów.")

    gmsh.initialize()
    gmsh.model.add("ground_rods")

    occ = gmsh.model.occ

    # Grunt:
    soil = occ.addBox(x_min, y_min, z_min, x_max - x_min, y_max - y_min, z_max - z_min )

    rod_volumes = []

    if len(rods) > 0:
        for i, ((d, l), (x, y)) in enumerate(zip(rods, xy), start=1):
            r = d / 2.0

            if z_max - l <= z_min:
                raise ValueError(f"Pręt {i}: długość l musi być mniejsza od grubości gruntu.")

            if x_min > x-r or x_max < x+r or y_min > y-r or y_max < y+r:
                raise ValueError(f"Pręt {i}: wychodzi poza obszar gruntu.")

            cyl = occ.addCylinder(x, y, z_max, 0.0, 0.0, -l, r)
            rod_volumes.append(cyl)

        cut, _ = occ.cut(
            [(3, soil)],
            [(3, tag) for tag in rod_volumes],
            removeObject=True,
            removeTool=True,
        )
        soil_volumes = [tag for dim, tag in cut if dim == 3]
    else:
        soil_volumes = [soil]

    occ.synchronize()

    gmsh.model.addPhysicalGroup(3, soil_volumes, 1)
    gmsh.model.setPhysicalName(3, 1, "soil")

    surfaces = gmsh.model.getEntities(2)

    eps = 1e-6

    top = []
    bottom = []
    x_min_surfaces = []
    x_max_surfaces = []
    y_min_surfaces = []
    y_max_surfaces = []

    for dim, tag in surfaces:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)

        if abs(zmin - z_max) < eps and abs(zmax - z_max) < eps:
            top.append(tag)

        elif abs(zmin - z_min) < eps and abs(zmax - z_min) < eps:
            bottom.append(tag)

        elif abs(xmin - x_min) < eps and abs(xmax - x_min) < eps:
            x_min_surfaces.append(tag)

        elif abs(xmin - x_max) < eps and abs(xmax - x_max) < eps:
            x_max_surfaces.append(tag)

        elif abs(ymin - y_min) < eps and abs(ymax - y_min) < eps:
            y_min_surfaces.append(tag)

        elif abs(ymin - y_max) < eps and abs(ymax - y_max) < eps:
            y_max_surfaces.append(tag)

    add_physical_surface("top", top, 101)
    add_physical_surface("bottom", bottom, 102)
    add_physical_surface("x_min", x_min_surfaces, 103)
    add_physical_surface("x_max", x_max_surfaces, 104)
    add_physical_surface("y_min", y_min_surfaces, 105)
    add_physical_surface("y_max", y_max_surfaces, 106)

    if len(rods) > 0:
        # Powierzchnie otworów po prętach
        for i, ((d, l), (x, y)) in enumerate(zip(rods, xy), start=1):
            r = d / 2.0
            rod_surfaces = []

            for dim, tag in surfaces:
                xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)

                inside_x = xmin >= x - r - 1e-5 and xmax <= x + r + 1e-5
                inside_y = ymin >= y - r - 1e-5 and ymax <= y + r + 1e-5
                inside_z = zmin >= -l - 1e-5 and zmax <= 1e-5

                if inside_x and inside_y and inside_z:
                    rod_surfaces.append(tag)

            add_physical_surface(f"rod_{i}", rod_surfaces, 200 + i)

        # Zagęszczanie siatki wokół prętów
        rod_surface_tags = []
        for i in range(1, len(rods) + 1):
            group = gmsh.model.getEntitiesForPhysicalGroup(2, 200 + i)
            rod_surface_tags.extend(group)
    else:
        rod_surface_tags = None

    if rod_surface_tags:
        gmsh.model.mesh.field.add("Distance", 1)
        gmsh.model.mesh.field.setNumbers(1, "FacesList", rod_surface_tags)

        gmsh.model.mesh.field.add("Threshold", 2)
        gmsh.model.mesh.field.setNumber(2, "InField", 1)
        gmsh.model.mesh.field.setNumber(2, "SizeMin", lc_rods)
        gmsh.model.mesh.field.setNumber(2, "SizeMax", lc_ground)
        gmsh.model.mesh.field.setNumber(2, "DistMin", 0.2)
        gmsh.model.mesh.field.setNumber(2, "DistMax", 2.0)
        gmsh.model.mesh.field.setAsBackgroundMesh(2)

    gmsh.option.setNumber("Mesh.MshFileVersion", msh_version)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", lc_rods)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", lc_ground)

    gmsh.model.mesh.generate(3)
    gmsh.write(output)

    gmsh.finalize()


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--a", type=float, default=None, help="szerokość X")
    parser.add_argument("--b", type=float, default=None, help="długość Y")
    parser.add_argument("--c", type=float, default=None, help="głębokość Z")

    parser.add_argument(
        "--box",
        type=parse_box,
        default=None,
        help='Zakres obszaru: "[x_min:x_max]x[y_min:y_max]x[z_min:z_max]"',
    )

    parser.add_argument(
        "--rods",
        type=str,
        required=False,
        help='Lista prętów w formacie średnica:długość "d:l,d:l", np. "0.016:0.8,0.02:1.2"',
    )

    parser.add_argument(
        "--xy",
        type=str,
        default=None,
        help='Pozycje prętów w formacie "x:y,x:y"; jeśli brak, pozycje automatyczne',
    )

    parser.add_argument("--lc-ground", type=float, default=5.0, help="gęstość całej siatki (domyślnie 5)")
    parser.add_argument("--lc-rods", type=float, default=0.01, help="gęstość w pobliżu prętów domyślnie (0.01)")
    parser.add_argument("-o", "--output", type=str, default="ground_rods.msh", help="nazwa pliku wynikowego (domyślnie ground_rods.msh")

    args = parser.parse_args()

    has_abc = args.a is not None or args.b is not None or args.c is not None
    has_box = args.box is not None

    if has_abc and has_box:
        parser.error("Podaj albo --a --b --c, albo --box, ale nie oba warianty.")

    if not has_abc and not has_box:
        parser.error('Musisz podać albo --a --b --c, albo --box "[x_min:x_max]x[y_min:y_max]x[z_min:z_max]".')

    if has_abc:
        if args.a is None or args.b is None or args.c is None:
            parser.error("Przy wariancie abc musisz podać wszystkie: --a --b --c.")

        x_min = 0.0
        x_max = args.a
        y_min = 0.0
        y_max = args.b
        z_min = -args.c
        z_max = 0.0
    else:
        x_min, x_max, y_min, y_max, z_min, z_max = args.box

    rods = parse_rods(args.rods) if args.rods is not None else []
    xy = parse_xy(args.xy) if args.xy is not None else None

    generate_mesh(
        x_min=x_min,
        x_max=x_max,
        y_min=y_min,
        y_max=y_max,
        z_min=z_min,
        z_max=z_max,
        rods=rods,
        xy=xy,
        output=args.output,
        lc_ground=args.lc_ground,
        lc_rods=args.lc_rods,
    )


if __name__ == "__main__":
    main()
