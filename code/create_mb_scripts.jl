using Pkg
Pkg.activate(".")
Pkg.instantiate()

##

datasets = readdir(
    "../data/correspondences_nexus",
    join=true
)

mkpath("mrbayes/output")

for nexus_file in datasets
    nm = split(split(nexus_file, "/")[end], ".")[1]*"_correspondences"
    mb = """
    #Nexus
    Begin MrBayes;
        execute ../$nexus_file;
        prset brlenspr = clock:uniform;
        lset rates=gamma;
        prset brlenspr = clock:uniform;
        prset clockvarpr=igr;
        lset coding=all;
        mcmcp stoprule=yes burninfrac=0.25 stopval=0.01 filename=output/$(nm) samplefreq=1000 printfreq=1000;
        mcmc ngen=100000000 nchains=4 nruns=4;
        sumt;
        sump;
        q;
    end;
    """
    write("mrbayes/$nm.mb.nex", mb)
end

##


datasets = readdir(
    "../data/cognate_classes_nexus",
    join=true
)



for nexus_file in datasets
    nm = split(split(nexus_file, "/")[end], ".")[1]*"_cognate_classes"
    mb = """
    #Nexus
    Begin MrBayes;
        execute ../$nexus_file;
        prset brlenspr = clock:uniform;
        lset rates=gamma;
        prset brlenspr = clock:uniform;
        prset clockvarpr=igr;
        lset coding=all;
        mcmcp stoprule=yes burninfrac=0.25 stopval=0.01 filename=output/$(nm)_cc samplefreq=1000 printfreq=1000;
        mcmc ngen=100000000 nchains=4 nruns=4;
        sumt;
        sump;
        q;
    end;
    """
    write("mrbayes/$nm.mb.nex", mb)
end
