#!/bin/sh

# bodge because the hudson_ms library stuff is not actually installed ...

# and to be able to use the mct library

# and to be able to use mrs library



#Edit this to match your path for the for hudson_ms_lib within your
#working copy of hudson_ms.

# 
#Edit the location and uncomment the line, and save this script as some other .sh filename
# (eg maybe 'custom_conf.sh')
# and use this 'personal' version of the script
#
# (you may have to change permissions so that you can execute the script:
# do this with chmod 755 custom_conf.sh 
# with your script file name for if it is not 'custom_conf.sh',
# or use a gui file manager to change the permissions.
#
# Doing this will ensure that the svn'd version of custom_config.sh 
# is not changed separately by everyone!

# ********** Set me! *************
#MS_DIR=/scratch-network/jah217/svn/trunk/hudson_ms/hudson_ms_lib

#MCT_DIR=/scratch-network/jah217/svn/trunk/mct/src

#MRS_DIR=/scratch-network/jah217/newsvn/mrs/branches/jah/src

#CXSCDIR=/usr/local/cxsc
#cxsc is included as isystem or we get too many warnings we can't do anything about...
set -x;

./configure CPPFLAGS="-I${MS_DIR} -I${MCT_DIR} -I${MCT_DIR}/aabc -I${MCT_DIR}/mct -I${MRS_DIR} -isystem${CXSCDIR}/include" LDFLAGS="-L${MS_DIR}/.libs -L${MCT_DIR} -L${MRS_DIR} -L${CXSCDIR}/lib" $@

