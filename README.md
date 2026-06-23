# lymph #
discontinuous po**LY**topal methods for **M**ulti-**PH**ysics &emsp; [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=lymphlib/lymph)

This repository contains a common interface for MATLAB packages using PolyDG methods with several application.

The full documentation of **lymph** is at [https://lymph.bitbucket.io/](https://lymph.bitbucket.io/).

**lymph** is developed at [MOX](http://mox.polimi.it) (The Laboratory for Modeling and Scientific Computing, Department of Mathematics) at Politecnico di Milano.
The development of **lymph** is part of the activities of the project **NEMESIS** ([https://erc-nemesis.eu](https://erc-nemesis.eu"), P.I.@PoliMI Paola F. Antonietti) and involves the following research group:

- Paola F. Antonietti (P.I.) <paola.antonietti@polimi.it>
- Stefano Bonetti <stefano.bonetti@polimi.it>
- Michele Botti <michele.botti@polimi.it>
- Mattia Corti <mattia.corti@polimi.it>
- Ivan Fumagalli <ivan.fumagalli@polimi.it>
- Ilario Mazzieri <ilario.mazzieri@polimi.it>

### Maintainers ###

* Stefano Bonetti <stefano.bonetti@polimi.it>
* Michele Botti <michele.botti@polimi.it>
* Mattia Corti <mattia.corti@polimi.it>
* Ivan Fumagalli <ivan.fumagalli@polimi.it>
* Ilario Mazzieri <ilario.mazzieri@polimi.it>

### How to include lymph in your code ###

If your MATLAB code (``OUTER_REPO`` from now on) is tracked by a git repository, the easiest way to include **lymph** in it (and to keep track of its updates) is by creating a git submodule:

0. If you have not already done it, set up a personal SSH key on your Bitbucket account (instructions [here](https://support.atlassian.com/bitbucket-cloud/docs/set-up-personal-ssh-keys-on-windows/); in most cases, you may start directly from point 3 of such instructions).
NB. If you have a Windows system, when generating the keys by `ssh-keygen` you may need to use `.ssh/id_rsa` as `{ssh-key-name}`, that is
    ```bash
    $ ssh-keygen -C "{username@emaildomain.com}" -f 'HOMEPATH/.ssh/id_rsa'
    ```
    where `HOMEPATH` is the full path to your home directory (e.g. `c/Users/yourusername`). See [here](https://stackoverflow.com/questions/20226147/where-does-github-for-windows-keep-its-ssh-key) for further info.
1. From the root directory ``OUTER_REPO`` of your code (in a clean state, i.e. with no pending modifications to be committed), add **lymph** as a submodule, by
    ```bash
    $ git submodule add git@bitbucket.org:lymph/lymph.git
    ```
    This will create a sub-directory ``lymph`` and a file ``.gitmodules``.
2. Commit the creation of the abovementioned sub-directory and file:
    ```bash
    $ git add .gitmodules lymph/
    $ git commit -m 'Add lymph to my repo'
    $ git push
    ```
3. In the main of your code, the following command adds all **lymph** code to your path:
    ```MATLAB
    addpath(genpath('lymph'));
    ```
    If you want to include only one module (e.g., `Core/MeshGeneration`):
    ```MATLAB
    addpath(genpath('lymph/Core/MeshGeneration'));
    ```
4. If you want to keep up-to-date with the last version of **lymph**, from time to time you need to update the submodule. From the root directory ``OUTER_REPO``:
    ```bash
    $ cd lymph
    $ git checkout main
    $ git pull
    $ cd ..
    $ git add lymph
    $ git commit -m 'Update lymph submodule'
    ```
5. If your code is installed in more than one location and you want to retrieve the same version of the submodule that you have elsewhere (e.g., after the update according to point 4), you can simply execute:
    ```bash
    $ git submodule update
    ```
    ***NB***: If you are updating the submodule for the first time in a local repository (e.g. you are installing your code on a new machine), the first time you need the ``--init`` option:

    ```bash
    $ git submodule update --init
    ```

The submodule is itself a git repository: in ``OUTER_REPO/lymph`` you can code, commit, create branches, and contribute to **lymph**!

### License ###

**lymph** is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License (see [LICENSE.md](LICENSE.md)), or (at your option) any later version. **lymph** is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. For more details, see the GNU Lesser General Public License included in this repository as [LICENSE.md](LICENSE.md) or visit the website [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/).

### Contribution guidelines ###

* `function`s should be named using the *CamelCase* convention.
* The blocks corresponding to different physics are ordered as: poroelasticity ('PE'), acoustics ('AC'), elasticity ('EL'), fluid ('FL'), prion spreading ('PS'), thermo-poroelasticity ('TPE'), multi-compartment poroelasticity ('MPE').
* The code documentation employs [Doxygen](https://www.doxygen.nl/), thanks to a MATLAB-C++ converter `m2cpp.pl` written in Perl (see [here](https://it.mathworks.com/matlabcentral/fileexchange/25925-using-doxygen-with-matlab)).
* To contribute with the addition of a new physics/problem, add a new subdirectory inside `Physics`.

### Useful links and tutorials ###

* `bash` command line interface on Unix/Linux ([tutorial](https://ryanstutorials.net/linuxtutorial/), [cheat sheet](https://github.com/RehanSaeed/Bash-Cheat-Sheet))
* [Markdown](https://bitbucket.org/tutorials/markdowndemo)
* [git](https://git-scm.com/) ([Pro Git book](https://git-scm.com/book/en/v2), [cheat sheet](https://training.github.com/downloads/github-git-cheat-sheet/), [visual cheat sheet](https://ndpsoftware.com/git-cheatsheet.html))
* [git submodules](https://www.atlassian.com/git/tutorials/git-submodule)
* Text editors: [Vim](https://www.vim.org/) (often opened by default at ``git commit``; [cheat sheet](https://devhints.io/vim)), [nano](https://www.nano-editor.org/) ([cheat sheet](https://www.nano-editor.org/dist/latest/cheatsheet.html))


### Acknowledgments ###

The developers of **lymph** acknowledge Prof. Paul Houston and Dr. Giorgio Pennesi for the original implementation of the quadrature-free approach adopted in the library [Antonietti, Houston, Pennesi. Journal of Scientific Computing 77 (2018)].

They acknowledge the support of

* ERC Synergy Grant n. 101115663 [NEMESIS: NEw generation MEthods for numerical SImulationS](https://erc-nemesis.eu)
* [ICSC-Centro Nazionale di Ricerca in High Performance Computing, Big Data, and Quantum Computing](https://www.supercomputing-icsc.it/) funded by European Union-NextGenerationEU
