lr2rmats_dir = ./lr2rmats
IsoModule_dir = ./IsoModule
rMATSEM_dir  = ./rmats-EM


# dependencies
#lr2rmats  = lr2rMATS this is just a bunch of deps in conda
IsoModule = ISOModule
rMATSEM   = rMATSEM


all:	 $(IsoModule) $(rMATSEM)

${IsoModule}:
	cd ${IsoModule_dir} && make || exit 255

${rMATSEM}:
	cd ${rMATSEM_dir} && Rscript install.R || exit 255

clean:
	rm -f ${lr2rmats_dir}/bin/* ${IsoModule_dir}/IsoModule

