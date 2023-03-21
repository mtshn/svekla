Software for predicting gas chromatographic retention indices and mass spectra.\
**Download here:** https://github.com/mtshn/svekla/releases  \
**Windows:** download and unpack release from this site, run **install.bat** (it will install Microsoft Visual C++ redistributable packages), then run **svekla.bat**\
Probably (if you have installed many other software before) the **svekla.bat** will work even without installation. Java is bundled to release!\
**Linux:** download and unpack release, run "svekla.sh". Writing to disk should be permitted. Java must be installed
(sudo apt-get install default-jdk)


CFM-ID, version 2.4 (https://sourceforge.net/projects/cfm-id/) is used for mass spectra. Models published at the links below are used for retention indices: \
https://doi.org/10.6084/m9.figshare.14602317 \
https://doi.org/10.6084/m9.figshare.12651680 \
\
This project is based on the code posted at these links.

The instructions how there releases were built are contained in the instructions_how_the_released_were_built.txt file. The procedure of compilation itself is very easy (install JDK, MAVEN and run "mvn package" command), but a lot of other binary dependencies and model files are required. The instructions_how_the_released_were_built.txt file explains were all these file were given.

The releases include many binary dependencies:

1) RDKit https://github.com/rdkit/rdkit  (3-clause BSD license)
2) lpsolve https://sourceforge.net/projects/lpsolve/ (LGPLv2 license)
3) JSME https://jsme-editor.github.io/ (3-clause BSD license)
4) CFM-ID 2.4  https://sourceforge.net/projects/cfm-id/  (LGPLv2 license)
5) Source code and pre-trained models from the above-mentioned links (MIT license)
6) Visual C++ Redistributable packages from Microsoft  (only Windows release)
7) Open JDK (only Windows release)
8) Deeplearning4j, JavaFX, and many other components with permissive and open source licenses were included to JAR file via Maven. See "pom.xml" for more information.


