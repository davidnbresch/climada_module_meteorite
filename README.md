climada_module_meteorites
=========================

climada module meteorites - a global meteorite hazard event set for climada

this repository contains an additional climada module, please install climada (the core module) first 
(see repository named simply climada)

in order to grant core climada access to additional modules, create a folder 'climada_additional' in the same directory as the core climada folder and copy/move any additional modules into climada_additional, without 'climada_module_' in the filename. E.g. if the addition module is named climada_module_MODULE_NAME, we should have
.../climada the core climada, with sub-folders as
.../climada/code
.../climada/data
.../climada/docs
.../climada_additional/MODULE_NAME with contents such as 
.../climada_additional/MODULE_NAME/code
.../climada_additional/MODULE_NAME/data
.../climada_additional/MODULE_NAME/docs
this way, climada sources all modules' code upon startup

see climada/docs/climada_manual.pdf to get started

copyright (c) 2014, David N. Bresch, david.bresch@gmail.com all rights reserved.
