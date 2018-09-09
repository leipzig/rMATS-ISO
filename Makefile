lr2rmats_dir = ./lr2rmats
IsoModule_dir = ./IsoModule
rMATSEM_dir  = ./rmats-EM


# dependencies
lr2rmats  = lr2rMATS
IsoModule = ISOModule
rMATSEM   = rMATSEM


all:		$(lr2rmats) $(IsoModule) $(rMATSEM)

${lr2rmats}:
	cd ${lr2rmats_dir} && make || exit 255

${IsoModule}:
	cd ${IsoModule_dir} && make || exit 255

${rMATSEM}:
	cd ${rMATSEM_dir} && Rscript install.R || exit 255

clean:
	rm -f ${lr2rmats_dir}/bin/* ${IsoModule_dir}/IsoModule

