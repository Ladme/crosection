# crosection: Cross-sectional area calculator

Calculates a cross-sectional area of a specified group of atoms in target dimension. Requires simulation frame or trajectory in which the selected group of particles does not cross the periodic boundaries!

## Dependencies

`crosection` requires you to have groan library installed. You can get groan from [here](https://github.com/Ladme/groan). See also the [installation instructions](https://github.com/Ladme/groan#installing) for groan.

## Installation

1) Run `make groan=PATH_TO_GROAN` to create a binary file `crosection` that you can place wherever you want. `PATH_TO_GROAN` is a path to the directory containing groan library (containing `groan.h` and `libgroan.a`).
2) (Optional) Run `make install` to copy the the binary file `crosection` into `${HOME}/.local/bin`.

## Options

```
Usage: crosection -c GRO_FILE -s SELECTION [OPTION]...

OPTIONS
-h               print this message and exit
-c STRING        gro file to read
-f STRING        xtc file to read (optional)
-n STRING        ndx file to read (optional, default: index.ndx)
-o STRING        output file (default: area.xvg)
-s STRING        selection of atoms which cross-section shall be calculated
-r STRING        file containing van der Waals radii of relevant atoms (default: vdwradii.dat)
-i CHAR          axis in which the cross-section shall be calculated (default: z)
-x FLOAT-FLOAT   grid dimension in x axis (default: box size from gro file)
-y FLOAT-FLOAT   grid dimension in y axis (default: box size from gro file)
-z FLOAT-FLOAT   grid dimension in z axis (default: box size from gro file)
-d INTEGER       density of the grid used for calculation in points per nm (default: 10)
```

Use [groan selection language](https://github.com/Ladme/groan#groan-selection-language) to select atoms. 

The `vdwradii` file (flag `-r`) must have the following format:
```
RESIDUE_NAME ATOM_NAME VDW_RADIUS
RESIDUE_NAME ATOM_NAME VDW_RADIUS
...
```
Blank lines in the file are ignored. See the provided `vdwradii_martini3.dat` for an example of a `vdwradii` file (containing Van der Waals radii for Martini 3 amino acids). Be extremely cautious when using the `vdwradii_martini3.dat` as terminal backbone beads of amino acids may have different sizes than usual! You may have to modify the `vdwradii_martini3.dat` according to your needs.

## Example usage

```
crosection -c md.gro -f md_centered.xtc -s Protein -n my_index.ndx -i x`
```

`crosection` will analyze the trajectory in `md_centered.xtc` obtaining the information about the character of atoms from `md.gro` (coordinates from `md.gro` are ignored). Average cross-sectional area of ndx group `Protein` will be calculated. The group `Protein` will be read from `my_index.ndx`. Van der Waals radii of the selected atoms will be read from `vdwradii.dat` (default option of the flag `-r`). The cross-sectional area of `Protein` will be calculated along the x-axis (flag `-i`) and written into an output file `area.xvg` (default option of the flag `-o`). `area.xvg` can be visualized using `xmgrace`.

## Limitations

DOES NOT TREAT PERIODIC BOUNDARY CONDITIONS PROPERLY! USE CENTERED SIMULATION TRAJECTORY.

Only tested on Linux. Probably will not work on anything that is not UNIX-like.
