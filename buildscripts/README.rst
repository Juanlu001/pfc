Build scripts
-------------

These are conda recipes needed to build the requirements for the project.

* AQUAgpusph https://github.com/sanguinariojoe/aquagpusph
* CBC.Solve https://github.com/Juanlu001/CBC.Solve (GitHub mirror)

To build the corresponding conda packages type::

    conda build aquagpusph cbcsolve --python 27
