### 1. Project Overview

This repository contains a specialized C++ package for algebraic graph theory, primarily built around the [Flag Algebra method] to compute bounds in extremal combinatorics.

This particular implementation focuses on using the core data structures and utility classes adapted from the original **Flag Algebra Program Package** to facilitate new research or applications in this area, particularly using local cuts computation, as introduced by Hu, Lidický, Martins and Volec in _Large multipartite subgraphs in h-free graphs_ (DOI: https://doi.org/10.1007/978-3-030-83823-2_113) and particularly on the approach outlined by Balogh, Clemen and Lidick-y in _Max Cuts in Triangle-free Graphs_ (DOI: https://doi.org/10.1007/978-3-030-83823-2_82). More details about it can be seen in my Undergraduate Thesis, available here, in Portuguese: .


### 2. Core Components Used

The following fundamental components were adapted and integrated from the original package:

* **Flag (`Flag.h`, `FlagTemplates.cpp`):** The main class representing a flag in the theory, including support for isomorphisms, normalization, and handling object types.
* **FlagGraph (`FlagGraph.h`, `FlagGraph.cpp`):** The underlying Universal Theory Object (UTO) used to represent graphs (or digraphs/tournaments, depending on the implementation) in the flag algebra framework. The I/O format for graphs from the original implementation is also used.
* **BigNum (`bignum.h`, `bignum.cpp`):** The arbitrary-precision arithmetic library necessary for exact rational calculations.

***

Most of additions are in the `tcc/` directory.

### 3. Acknowledgment and Attribution

This work is based on the **Flag Algebra Program Package**, originally developed and copyrighted by **Leonardo Nagami Coregliano**.

The core components for flag algebra implementation, specifically the `Flag`, `FlagGraph`, and `BigNum` classes, are derived from this seminal project.

**Original Project Details:**

* **Link to GitHub repository:** https://github.com/robertoparente/flag-algebra-program-package.git.
* **Project Name:** Flag algebra semidefinite CPP package.
* **Primary Author & Copyright:** Leonardo Nagami Coregliano, Copyright (C) 2018.
* **Associated Research:** The original package was used to obtain results for the paper "On the maximum density of fixed strongly connected subtournaments," coauthored with Roberto Freitas Parente and Cristiane Maria Sato.

***

### 4. Licensing

This program is distributed under the terms of the **GNU General Public License, Version 3 (GPLv3)**.

The entire work is licensed under **GPLv3** to comply with the license requirements of the original source code from which the core components were derived.

You are free to redistribute and/or modify this program under the terms of the GPLv3, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but **WITHOUT ANY WARRANTY**; without even the implied warranty of **MERCHANTABILITY** or **FITNESS FOR A PARTICULAR PURPOSE**.

A copy of the GNU General Public License is included in the file `GPL.txt` and can be found online at [https://www.gnu.org/licenses/](https://www.gnu.org/licenses/).


