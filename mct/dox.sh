#!/bin/sh

#Make and run my own doxygen configuration file

#Bases configuration on current Doxygen file and adapts it
#Shows internal comments
#Takes out all the graphs stuff

DOXYFILE=Doxyfile
MYDOXYFILE=myDoxyfile
CMDFILE=sed_cmd_file

#replace as per the command file
sed -f ${CMDFILE} ${DOXYFILE} > ${MYDOXYFILE}

# run doxygen
doxygen ${MYDOXYFILE}

