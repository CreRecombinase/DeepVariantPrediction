#configfile: "config.json"
include: "helpers.py"

                

rule all:
    input:"ASD_SNV_peaks.h5"



rule subsample_h5:
    "This also binarizes the output"
    input:
        h5file = "deepsea_train/{project}_train.mat"
    params:
        nsamples=150000,
        newoutput=0
    output:
        h5file = "deepsea_train/{project}_neg_train.mat"
    run:

        hf = h5py.File(input.h5file,'r') # HDF5 file
        hf_labels = hf['/traindata']
        hf_data = hf['/trainxdata']
        hshape = hf_data.shape
        hf.close()
        nsamples = int(params.nsamples)
        print("Generating choice!")
        sub_sample=np.random.choice(hshape[2]-1,nsamples)
        sub_sample=np.sort(sub_sample)
        i=0
        for s_ind in grouper(sub_sample,20000):
            print(i)
            t_data=read_array_index(input.h5file,
                                    '',
                                    'trainxdata',
                                    np.array(s_ind))
            s_labels=np.zeros((1,t_data.shape[2]))
            write_array_h5(output.h5file,"/","trainxdata",t_data)
            write_mat_h5(output.h5file,'/',"traindata",s_labels,axis=1)
            i=i+1
               
    
rule center_bed:
    input:
        bedfile = "{project}.bed"
    params:
        bed_windowsize=100
    output:
        bedfile = "{project}_centered.bed"
    run:
        coord_f=pd.read_csv(input.bedfile,
                            sep="\t",
                            header=None,
                            names=['chrom','start','end'],
                            usecols=[0,1,2])
        middle=((coord_f['end']+coord_f['start'])/2).round(0).astype('int64')
        coord_f['start']=middle-int(params.bed_windowsize)
        coord_f['end']=middle+int(params.bed_windowsize)
        coord_f.ix[coord_f.start<0,'start']=0
        coord_f.to_csv(output.bedfile,sep="\t",header=None,index=False)
        
        

rule extend_bed:
    input:
       bedfile ="{project}_centered.bed"
    params:
        rev_comp=True,
        extend_length=400
    output:
        outbedfile="{project}_extend.bed"
    run:
        import pandas as pd

        coord_f=pd.read_csv(input.bedfile,
                            sep="\t",
                            header=None,
                            names=['chrom','start','end'],
                            usecols=[0,1,2])
        coord_f['start']=coord_f['start']-int(params.extend_length)
        coord_f['end']=coord_f['end']+int(params.extend_length)
        tmed=coord_f.ix[coord_f.start<0,'end']-coord_f.ix[coord_f.start<0,'start']
        coord_f.ix[coord_f.start<0,'end']=tmed
        coord_f.ix[coord_f.start<0,'start']=0
        coord_f['strand']='+'
        rcoord=coord_f
        rcoord['strand']='-'
        frcoord=pd.concat([coord_f,rcoord])
        frcoord.to_csv(output.outbedfile,sep="\t",header=None,index=False)

rule bed_to_fasta:
    input:
        bedfile= "{project}_extend.bed",
        hg19="../hg19/hg19.fa"
    output:
        fastafile="{project}_extend.fasta"
    shell:
        "bedtools getfasta -fi {input.hg19} -s -bed {input.bedfile} > {output.fastafile}"

rule train_test_validate:
    input:
        posfile="{project}_peaks.h5",
        negfile="deepsea_train/{project}_neg_train.mat"
    params:
        test_size=0.33,
        seed=1337
    output:
        trainfile="{project}_train.h5",
        testfile="{project}_test.h5",
        validfile="{project}_valid.h5"
    run:
        from sklearn.model_selection import train_test_split
        negmat = h5py.File(input.negfile,'r')
        posmat = h5py.File(input.posfile,'r')

        X_neg = np.transpose(np.array(negmat['trainxdata']),
                             axes=(2,0,1))
        X_pos =np.transpose(np.array(posmat['trainxdata']),
                            axes=(2,0,1))

        X_mat = np.concatenate([X_neg,X_pos],axis=0)
        
        y_neg= np.array(negmat['traindata'])
        y_pos = np.array(posmat['traindata'])
        
        y_mat = np.concatenate([y_neg,y_pos],axis=1).T
        
        X_train_m, X_test, y_train_m, y_test = train_test_split(X_mat,
                                                            y_mat,
                                                            test_size=0.15,
                                                            random_state=1337)
        X_train, X_valid, y_train, y_valid = train_test_split(X_train_m,
                                                            y_train_m,
                                                            test_size=0.15,
                                                              random_state=1337)
        
        
        write_array_h5(output.trainfile,"","trainxdata",X_train)
        write_mat_h5(output.trainfile,"","traindata",y_train)
        write_array_h5(output.testfile,"","testxdata",X_test)
        write_mat_h5(output.testfile,"","testdata",y_test)
        write_array_h5(output.validfile,"","validxdata",X_valid)
        write_mat_h5(output.validfile,"","validdata",y_valid)
        
        

        

rule fasta_to_h5:
    input:
        fastafile="{project}_extend.fasta"
    params:
        chunksize=10000
    output:
        h5file="{project}_peaks.h5"
    run:

        tseq = ['a','g','c','t','n']
        tord = [ord(el) for el in tseq]
        i=0
        ogrp='/'
        datan='trainxdata'
        print(tord)
        for chunk in grouper(fastest_fasta(input.fastafile),int(params.chunksize)):
            print(i)
            print("Converting")
            print(len(chunk))
            print(chunk[0][0:5])
            bin_mat = np.array(
                [
                    preprocessing.label_binarize(
                        [ord(el) for el in fasta_seq.lower()],
                        classes=tord[0:4])for fasta_seq in chunk if fasta_seq!=None ]
                )
            print(bin_mat.shape)
            print("Transposing")
            bin_mat=bin_mat.transpose((1,2,0))
            
            res_mat=np.ones((1,bin_mat.shape[2]))
            print("Writing")
            write_array_h5(output.h5file,"","trainxdata",bin_mat)
            write_mat_h5(output.h5file,"/","traindata",res_mat,1)
            print("Reading")
            i=i+1

        
