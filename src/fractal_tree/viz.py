from __future__ import annotations
import xml.etree.cElementTree as ET
from typing import Union
from pathlib import Path
import numpy as np


def write_line_VTU(
    nodes: np.ndarray, elements: list[tuple[int, int]], filename: Union[Path, str]
) -> None:
    file = ET.Element("VTKFile")
    file.set("type", "UnstructuredGrid")
    file.set("version", "0.1")
    file.set("byte_order", "BigEndian")
    UG = ET.SubElement(file, "UnstructuredGrid")
    piece = ET.SubElement(UG, "Piece")
    piece.set("NumberOfPoints", str(len(nodes)))
    piece.set("NumberOfCells", str(len(elements)))

    points = ET.SubElement(piece, "Points")
    DAp = ET.SubElement(points, "DataArray")
    DAp.set("type", "Float32")
    DAp.set("NumberOfComponents", "3")
    DAp.set("Format", "ascii")
    DAp.text = ""
    DAp.text = "\n".join(
        map(lambda a: str(a[0]) + " " + str(a[1]) + " " + str(a[2]), nodes)
    )

    cell = ET.SubElement(piece, "Cells")
    DAc = ET.SubElement(cell, "DataArray")
    DAc.set("type", "Int32")
    DAc.set("Name", "connectivity")
    DAc.set("Format", "ascii")
    DAc2 = ET.SubElement(cell, "DataArray")
    DAc2.set("Name", "types")
    DAc2.set("Format", "ascii")
    DAc2.set("type", "Int32")
    DAc3 = ET.SubElement(cell, "DataArray")
    DAc3.set("Name", "types")
    DAc3.set("type", "Int32")
    DAc3.set("Format", "ascii")
    DAc3.set("Name", "offsets")
    DAc.text = ""
    DAc2.text = ""
    DAc3.text = ""
    DAc.text = "\n".join(map(lambda a: str(a[0]) + " " + str(a[1]), elements))
    DAc2.text = "\n".join(["3"] * len(elements))
    DAc3.text = "\n".join(map(str, range(2, len(elements) * 2 + 1, 2)))
    tree = ET.ElementTree(file)
    tree.write(filename)
