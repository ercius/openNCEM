# ncempy.fio

The `ncempy.fio` module contains the file IO necessary for the various dataformats floating around in electron microscopy. Internally the package is designed to work in the versatile EMD file format, other formats are interfaced with importers. The interfaces of each file format are implemented in their own classes.

## Contents

Overview of contents with short description:

| Submodule | Description |
| --- | --- |
| `emd` | EMD file format |
| `ser` | SER file format used by TIA (FEI) |
