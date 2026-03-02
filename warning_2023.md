### *WARNING (31th March 2023)*
It seems that the various Jupyter Notebooks which use the ipympl library have a problem with the latest versions of the anaconda distribution (starting with `python 3.9`), and probably with other distribution as well.

Try the following possibilities if you have problem with the graphic and/or the interaction in the Spike Notebooks

1. the team developping jupyter is migrating the notebook to a new technology and Spike is not yet adapted to this new environment.
   In the meantime, use `nbclassic` instead of `notebook` to launch the various notebooks
    ```
    jupyter-nbclassic TheNooteBook.ipynb
    ```
2. The interaction is based on the optionnal `ipympl` extension, and some people have instabilities in the interaction
*upgrading ipympl to the 0.9.3 version seems to solve the difficulty:*
    ```
    pip install ipympl==0.9.3
    ```
or alternatively
    ```
    conda config --env --add channels conda-forge
    conda install ipympl=0.9.3
    ```
Depending on your set-up (jupyter version) it might work either in `jupyter lab`,  `jupyter notebook` or `jupyter nbclassic`, or even all of them.


