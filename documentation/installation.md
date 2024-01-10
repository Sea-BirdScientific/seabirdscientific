# Installing the Sea-Bird Scientific Community Toolkit

<!-- TOC here
[VSCode and Python Setup](VSCode and Python Setup) 
-->

## Using a Local Checkout of the Toolkit Repository and VSCode
### VSCode and Python Setup
If you have [VSCode](https://code.visualstudio.com/) and [Git](https://www.git-scm.com/downloads) already installed on your computer:
1. Checkout the [Sea-Bird Scientific](https://github.com/Sea-BirdScientific/seabirdscientific.git) repository from GitHub using a Git GUI Client or one of the following command line methods:
    - To place the toolkit in a folder named seabirdscientific:
        1. At a command prompt, change to the parent directory where you would like the toolkit to reside.
        1. Run the command ```git clone https://github.com/Sea-BirdScientific/seabirdscientific.git```
    - To place the toolkit in a folder with a different name:
        1. At a command prompt, change to the parent directory where you would like the toolkit to reside.
        1. Run the command ```git clone https://github.com/Sea-BirdScientific/seabirdscientific.git yourFolderName```
1. Run VSCode and choose File | Open Folder. 
1. Use the file picker to select the toolkit folder that was just created.
1. Use the [VS Code - Creating Environments](https://code.visualstudio.com/docs/python/environments#_creating-environments) page as a guide to setup a virtual environment for your work.
1. Open a VSCode power shell terminal if one is not currently active
1. Ensure that the virtual environment by hovering over the terminal tab and checking the Python: extension for the activated environment
<img src="images/VerifyEnv.PNG" width=400>
1. Execute the following command in the powershell terminal to install the toolkit and its dependencies
    ``` bash
    py -m pip install -e .[dev]
    ```
1. Open the file documentation/ctd-processing.ipynp notebook in VSCode and run the cells to see the notebook in action.
    
### VSCode and Anaconda/Miniconda

## Related Software
1. [Python](https://www.python.org/downloads/) version 3.9 or higher
1. [Git](https://www.git-scm.com/downloads) 
1. [Git GUI clients](https://www.git-scm.com/downloads/guis)
1. [Anaconda]()
1. [Visual Studio Code](https://code.visualstudio.com/)
    1. Install and start [VSCode](https://code.visualstudio.com/)
    1. Type Ctl-Shift-X to open the extension manager
       1. Ensure that the Python extension from Microsoft is installed 
       1. Ensure that the Jupyter extension from Microsoft is installed
