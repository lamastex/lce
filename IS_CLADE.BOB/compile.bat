cl /Feconstant.exe -DCONSTANT is.c uniform.c
cl /Feexponential.exe -DEXPONENTIAL is.c uniform.c
cl /Feexp_flat.exe -DEXP_FLAT is.c uniform.c
del *.obj
