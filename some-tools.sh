# seqtk - getting fastq sequences for each cluster
fseqtk(){
    cd ${BASEDIR}/${WORKDIR}
    for file in ${BASEDIR}/${WORKDIR}/02-400-do-assembly/*
    do
        for seq in $file
        do
            /home/tools/seqtk/seqtk subseq ${BASEDIR}/sequencing_data/barcode02.fastq $seq > ${BASEDIR}/${WORKDIR}/seqtk-09-400/$(basename "$seq").fastq
        done
    done
}

# Flye - assembly of each cluster
f-flye-klastrow(){
    cd ${BASEDIR}/${WORKDIR}
    for klaster in ${BASEDIR}/${WORKDIR}/seqtk-02-400/*
    do
        /home/tools/Flye-2.9.2/bin/flye --threads 40 --nano-hq $klaster --out-dir ${BASEDIR}/${WORKDIR}/flye-02-700-2/$(basename "$klaster")
    done
}

# Blat - alignment of a gene sequence against assembly sequences
fblat-assembly(){
	cd ${BASEDIR}/${WORKDIR}
	for klaster in ${BASEDIR}/${WORKDIR}/flye-04-400-2/*
	do
        # searching for assembly file in current directory
        CHECKING=$(find "$klaster" -maxdepth 2 -name "assembly.fasta")
        # alignment 
		/home/tools/blat/blat -noHead /data/collaborators/Phage_Antigen_Presentation/gene3.fasta ${CHECKING} ${klaster}/gene3.psl
	done
}