# ITG Challenge

## Description/Instructions

Please find below a few questions that mimic some common problems we encounter at ITG. They are grouped by broad theme. You will notice there are many questions, the goal is not to answer them all but to pick a few questions to focus on (10 is a good number, but pick as many as you want). You should pick from all three categories, but there are many more bioinformatics questions so you should mainly pick from those. We encourage you to choose your questions according to your areas of expertise but also to try and answer questions that are as varied as possible.

For programmatic questions, you can use the language and libraries of your choice, but we will assess whether your choice of language was optimal. Try and aim for a minimal solution in terms of code length. If you use a shell script, you can assume that common non-core packages will be installed (e.g. `awk`, `sed`, `perl`, `python`, `sponge`, `wget` or `jq`). You can use the shell of your choice, if not otherwise specified we will assume `bash`. Assume that all common bioinformatics tools `bcftools`, `bedtools`, `vcftools`, `plink` and others are all installed.

We are primarily interested in how you would solve these problems if you encountered them in real life. Whenever the command line or programming is used, please include your code along with your answer. Not all questions are programmatic, and some are open-ended. Feel free to include details and to discuss potential issues if you don't think there is a clear-cut answer.

To submit your results, please clone this repository and make your edits. Once you're done, send us a link to your repository, or compress it and send us a link to the archive.

## Questions

### Support/resource management/Shell
1. A user has several versions of R installed in their path. Each version of R has a number of locally installed libraries. The user is confused, and would like to know which library is installed for each version of R. Can you write a command to help them out?
    <ul>
    <b>Ans:</b> The path to R libraries should be present in .libPaths() and `installed.packages()` can be used to get the information about the installed packages in the default library path which is .libPaths()[1] by default. So, by running `lapply(.libPaths(), installed.packages)` the user can get package information from the different library paths.
    </ul>

3. A common problem with shared filesystems is a disk quota overflow. This can be due to 1) a large combined size on disk or 2) too many files present on disk. We would like to help users who encounter this problem to locate the problematic files/directories. Write a command to sort all subdirectories of level `n` (`n` determined by the user) by their human-readable size. Write another command to sort all subdirectories of level `n` according to the number of files they contain.
4. A user wants to install an `R` package and gets the following [error log](data/error.log). What is likely to cause the error and how can they solve it?
    <ul>
    <b>Ans:</b> The error log contains the following error several times:<code>cc1plus: error: unrecognised command line option ‘-std=c++11’</code>. This points to an incompatible g++ compiler, perhaps an older one which does not recognize the  `-std=c++11` flag. Also, towards the end the error like `/bin/sh: nc-config: command not found` suggests that some of the development package is missing or not in the PATH. So, the user needs to check these points, and fixing them might solve the issue.
    </ul>
6. A user is running commands like this one `cat file1 <(cut -d " " -f 1-15,17,18 file2) > file3`. What does this command do? It runs fine on the command line, but then the user includes it into a file with other commands, saves it and runs `chmod +x` on it. However, that line of code throws the following error : `syntax error near unexpected token '('`. What has the user forgotten?
    <ul>
    <b>Ans:</b> The command extracts the columns 1st to 15th, then 17th, 18th from file2 and then vertically concatenates at the end of file1 and finally saves into file3. I think the user forgot to include the proper shebang line at the top of the script. By adding a line like <i>#!/bin/bash</i> should solve the issue.
    </ul>
8. A collaborator has sent you [this script](data/EasyQCWrapper.sh). It is a wrapper for a bioinformatics software called `EasyQC`.  Running it, you get the following error: 

    ```bash
    ./test.EasyQC-START.R: line 6: syntax error near unexpected token 'EasyQC'
    ./test.EasyQC-START.R: line 6: 'library(EasyQC)'
    ```

     You need to run this script now, but your collaborator is unavailable for a few days. What is causing the error? (Hint: Nothing is wrong with the `.ecf` EasyQC script.)
9. Programmatic download
    - You have to download all autosomal files from this location: [http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/) onto **your server**. You connect to the server via SSH. Using only the command line, how do you perform this download?
    - You are at a conference abroad and you quickly realise that your connection is unstable. You get disconnected constantly, which interrupts the download. How do you ensure the download survives these disconnections?
10. Bioinformaticians often work on a computing cluster. The cluster runs a software called a job scheduler that attributes resources to users depending on the requirements of their jobs. In this case, let's imagine the cluster is running IBM LSF. You do not need to know it to answer this question. The `bjobs` command lists all jobs submitted by the user (manpage [here](https://www.ibm.com/support/knowledgecenter/en/SSETD4_9.1.2/lsf_command_ref/bjobs.1.html)). It gives this kind of output:
    ```
    JOBID   USER             STAT  QUEUE      FROM_HOST EXEC_HOST JOB_NAME SUBMIT_TIME
    9670681 current_user     RUN   basement   head_node node1     job1     Oct 24 10:24
    9740051 current_user     RUN   basement   head_node node1     job2     Oct 24 17:41
    9670681 current_user     RUN   normal     head_node node2     job3     Oct 24 10:24
    9740981 current_user     PEND  basement   head_node           job4     Oct 24 17:44

    ```
     - Given the [following output](data/farm-snapshot.txt) of `bjobs -all`, which users are the top 5 users of the cluster?
     - How many jobs does the user `pathpip` have running in all queues?
     - A user wants to know how many jobs they have pending (`PEND`) and running (`RUN`) in each queue. Write a command line to do that (You can use the log above to check your command line). How would they display this on their screen permanently, in real time?
11. An analysis you need to run on the cluster requires a particular python library, but you do not have administrator rights. IT is on holiday. What do you do?
    <ul>
    <b>Ans:</b> I will install it locally using `pip install --user my_package`. This will install it in ~/.local/lib/pythonX.Y/site-packages/my_package/, where X, and Y are major and minor Python versions. Then, I will also need to add the ~/.local/bin to PATH as export PATH=$HOME/.local/bin:$PATH in my .bashrc file. This is done to make sure my package binary is discoverable during runtime in case it is a command-line tool.
    </ul>

12. All major computational tasks in your lab are done via SSH connection to mainframe servers or HPC clusters. A user comes from a Linux (mostly command-line) background but IT only support Windows 10 for laptops. How would you advise them to configure their laptop to make their transition easier?
    <ul>
    <b>Ans:</b> In my opinion, with the following configurations, a user coming from a Linux background will find Windows 10 much more comfortable for interacting with HPC:
	    <ul>
		    <li>WSL Linux: By installing WSL, e.g., Ubuntu on Windows, the user can get a near‐native Linux environment and command line.</li>
		    <li>Git Bash: This also offers Linux like shell on Windows. However, installing utility tools like <i>rsync</i> is a bit involved.</li>
		    <li>VS Code: Optionally, if the user is used to VS code IDE, then using that on Windows can boost productivity as it can support both WSL and Git bash terminals.</li>
	    </ul>
    </ul>

### Bioinformatics
1. The [VCF format](http://www.internationalgenome.org/wiki/Analysis/vcf4.0/) is a popular format to describe genetic variations in a study group. It is often used in sequencing projects. Due to size concerns, it is often compressed using `gzip` and indexed using `tabix`. A binary version, BCF, also exists.
    - Write a command or script to remove duplicate positions in a VCF such as [this one](data/duplicates.vcf.gz), independently of their alleles. The positions can be duplicated an arbitrary number of times. Write code to keep the first, last and a random record among each set of duplicated records.
    - Same question, but make duplicate detection allele-specific. When it finds such an exact duplicate, your code should remove all of the corresponding records.
    <p><b>Ans:</b></p>
      
    ```
    import pysam
    import random
    import pdb
    import pandas as pd
    import gzip
    
    ## User defined functions
    def select_random_without_first_last(group):
        """Selects a random row excluding the first and last records."""
        if len(group) <= 2:
            return None  # No middle records available
        middle_records = group.iloc[1:-1]  # Exclude first and last
        return middle_records.sample(1)
    
    # Read the VCF file and extract variant records while handling missing columns
    vcf_file_path = "duplicates.vcf.gz"
    variant_records = []
    header_lines = []
    
    with gzip.open(vcf_file_path, "rt") as vcf_file:
        for line in vcf_file:
            if line.startswith("#"):
                header_lines.append(line)  # Store header for later writing
            else:
                fields = line.strip().split("\t")
                variant_records.append(fields)
    
    # Determine the maximum number of columns present
    max_cols = max(len(record) for record in variant_records)
    
    # Define flexible column names based on VCF specification
    vcf_columns = header_lines[-1].strip("\n").replace("#", "").split("\t")
    
    # Convert to DataFrame
    df = pd.DataFrame(variant_records, columns=vcf_columns[:max_cols])
    
    # Convert numeric columns to appropriate types
    df["POS"] = df["POS"].astype(int)
    
    # Task 1. keep first, last, and a random entry
    df_sorted = df.sort_values(by=["CHROM", "POS"])
    df_grouped = df_sorted.groupby(["CHROM", "POS"])
    # pdb.set_trace()
    
    df_first = df_grouped.first().reset_index()
    df_last = df_grouped.last().reset_index()
    # df_random = df_grouped.apply(lambda x: x.sample(1)).reset_index(drop=True)
    df_random = df_grouped.apply(select_random_without_first_last).reset_index(drop=True)
    
    df_position_filtered = pd.concat([df_first, df_last, df_random]).drop_duplicates()
    
    # 2. Allele-specific duplicate removal (Remove exact duplicates)
    df_allele_grouped = df_sorted.groupby(["CHROM", "POS", "REF", "ALT"])
    df_allele_filtered = df_allele_grouped.filter(lambda x: len(x) == 1)  # Remove exact duplicates
    
    # Save both results as new VCF files
    position_filtered_path = "cleaned_positions.vcf.gz"
    allele_filtered_path = "cleaned_allele_specific.vcf.gz"
    
    # Writing the results back to VCF format
    def write_vcf(output_path, df_filtered):
        with gzip.open(output_path, "wt") as vcf_out:
            vcf_out.writelines(header_lines)  # Write header first
            for _, row in df_filtered.iterrows():
                vcf_out.write("\t".join(map(str, row.tolist())) + "\n")
    
    # Save cleaned VCFs
    write_vcf(position_filtered_path, df_position_filtered)
    write_vcf(allele_filtered_path, df_allele_filtered)
    ```
2. From an existing VCF with an arbitrary number of samples, how do you produce a VCF file without any samples using `bcftools`?
   <p><b>Ans:</b></p>
   
    ```
    bcftools view -G samples.vcf.gz -Oz -o no_samples_my.vcf.gz
    ```
4. You are the curator of a genotype dataset with a very strict privacy policy in place. In particular, it should be impossible to tell, given access to a person's genetic data, whether they were part of your study by looking at a dataset you provided. A collaborator is asking you for some data to run tests on their code. What information can you safely contribute from your study?
5. How do you convert a gzipped VCF to the `bimbam` format? (you may choose to script a solution yourself, or not)
6. A user sends you a small number of chromosome and positions in build 38 that they want to know the rsID of. 
    - What is missing from their request? What kind of unexpected output can they expect?
    - Given [this file](data/rand.chrpos.txt), honour their request using the Ensembl REST API.
    - Do the same, but offline, using the dbSNP r.150 VCF file.
    - What would change if these positions were in build 37?
    - If the user sends you 7,000 such chromosome positions, how would the above methods perform? Do you know of any alternatives?
7. How would you change the chromosome numbers in the file above to chromosome names (e.g. "chr1" instead of "1")?
    - How would you change the names back to the original? Would your solution work if an additional column containing text of arbitrary length and content is appended at the left of the file?
    - These positions are extracted from a VCF. Convert this file to the BED format.
      <p><b>Ans:</b> To change the chromosome names (e.g. "chr1" instead of "1"),  I will do <code>awk '{ $1 = "chr"$1 }1' rand.chrpos.txt > rand.chrpos_with_chr.txt</code>. However, if a new column is inserted then the previous solution may not work depending on where the column is inserted. In that case, we need to adjust the solution to match where <i>chromosome</i> and <i>position</i> appear. For example, if a column is added on the left then the modified snippet would be <code>awk '{ $2 = "chr"$2 }1' rand.chrpos.txt > rand.chrpos_with_chr.txt</code> and similarly to revert it to the original format we can use <code>awk '{ sub(/^chr/,"",$2) }1' rand.chrpos_with_chr.txt > rand.chrpos_without_chr.txt </code>. Finally, to convert to BED which is 0-based, half-open intervals, I can do as <code>awk 'BEGIN { print "chrN", "start", "end"} { print "chr"$1, $2-1, $2 }' OFS="\t" rand.chrpos.txt > rand.chrpos.bed</code>
      </p>
8.	Download the 1000 Genomes sites VCF file for chromosome 21 [here](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr21_GRCh38_sites.20170504.vcf.gz). We want to compare it to [a locally stored file](data/ALL.chr21_GRCh38_sites.20170504.vcf.gz).
    - What is the fastest way to check the integrity of, or compare, any such downloaded file?
    - If you find that the files are indeed different, how do you find their differences? Keep in mind that this kind of file can be very large (>100Gb), your solution should be as fast and memory-efficient as possible.
    - If you found no differences in the end, what could cause a false alarm?
    <p></p>
    <ul>
    <b>Ans:</b>
	<ul>
	<li> The given link does not work so I will provide a conceptual answer here. I will use <i>md5sum</i> to compute the checksums of the files and compare them using <i>diff</i> command. If no output is produced then the files are identical.</li>
    <li>To compare large VCF files, one can use <i>bcftools</i>, e.g.,<code>bcftools cmp remote/ALL.chr21_GRCh38_sites.20170504.vcf.gz local/ALL.chr21_GRCh38_sites.20170504.vcf.gz</code>. This generates a summary of matching vs. non-matching records and is typically more memory‐efficient than a raw line‐by‐line <i>diff</i>, because <i>bcftools</i> can parse BGZF blocks and parse VCF records in a streamed fashion.</li>
    <li>In the event of false alarm, you may find the compressed files differ in a bitwise sense but the actual VCF contents might still be functionally identical. Some common reasons could be, different compression, or minor changes in header lines (dates, software versions), or other metadata differences that don’t affect the actual genotype data.</li>
    </ul>
</ul>
9.	What is the p-value corresponding to standard normal z-scores of 10.35, 29.7, 45.688 and 78.1479?
    <ul>
    <b>Ans:</b> Since all given Z-scores are positive, we only need to compute the right-tailed p-values (i.e., 𝑃(𝑍 > 𝑞)).
    In R this can be done by using the following commands,
    <br>q <- c(10.35, 29.7, 45.688, 78.1479) <br>
    p_values <- pnorm(q, lower.tail = FALSE) # One-tailed (right-sided) p-values
    </ul>    

11.	We want to round a column of numbers to `n` decimal places, with values with 5 as their rightmost significant digit rounded up. Use the language of your choice.
     <p><b>Ans:</b></p>
     
    ```
    import pandas as pd
    from decimal import Decimal, ROUND_HALF_UP
    import pdb
    
    def round_half_up(x, n=0):
        """
        Rounds a float x to n decimals. If the rightmost digit in the
        (original) decimal part is '5', then the nth digit is bumped up 
        (with proper carry propagation) regardless of any conventional
        rounding rules.
    
        Parameters:
          x: the number (typically a float) to round.
          n: the number of decimal places desired.
        
        Returns:
          The rounded value as a float.
        """
        # Work only if x is a float (otherwise, leave it unchanged)
        if isinstance(x, float):
            s = str(x)  # use the float’s string representation
            if '.' not in s:
                return x  # no fractional part
            
            integer_part, frac_part = s.split('.')
            
            # If the caller asked for more decimals than exist, use what we have.
            if n > len(frac_part):
                n = len(frac_part)
            
            # If there are no decimals requested, work on the integer part.
            if n == 0:
                # Check if the entire number ends with a 5.
                if s[-1] == '5':
                    return float(int(integer_part) + 1)
                else:
                    return float(round(x))
            
            # If the original fractional part (as given) ends with '5'
            # then we force an upward adjustment in the target n-digit portion.
            if frac_part[-1] == '5':
                # pdb.set_trace()
                # Work with the first n digits (even if frac_part has extra digits,
                # we only consider the desired precision)
                target = list(frac_part[:n])
                # Increment the last digit of the target.
                # (This is the “round up” step.)
                # We must handle the possibility of a carry if the digit is 9.
                i = n - 1
                carry = 1
                while i >= 0 and carry:
                    new_val = int(target[i]) + carry
                    if new_val == 10:
                        target[i] = '0'
                        carry = 1
                    else:
                        target[i] = str(new_val)
                        carry = 0
                    i -= 1
                # If a carry remains, increment the integer part.
                if carry:
                    integer_part = str(int(integer_part) + 1)
                new_frac = "".join(target)
                return float(integer_part + '.' + new_frac)
            else:
                # Otherwise, use a conventional round()
                return round(x, n)
        else:
            return x
    
    # Example usage
    if __name__ == '__main__':
        # n = 14
        n = 2
        pd.set_option("display.precision", 14)
        df = pd.DataFrame({
            'values': [2.41238537658565, 1.225658565, 1.230, 1.2005, 1.225]
        })
        df['rounded'] = df['values'].apply(lambda x: round_half_up(x, n))
    ```
13.  Is [this HRC-imputed file](https://drive.google.com/open?id=1dOYlwIlAyz9-i4gVy2sgpQv_4YX9Y5nU) missing any chromosomes? Try to find out in seconds if you can.
14.  Find out the coordinates of the _ADIPOQ_ gene. Your method should be generalisable to a list of genes and take a few seconds to run (take inspiration from question 5). Then, please report the following:
    - the coordinates of the exons of its canonical transcript.
    - all documented variants in this gene.
    - all phenotype-associated variants. 
    - all documented loss-of-function (LoF) variants in this gene. How did you define LoF?
    - If asked to find all regulatory variants that potentially affect this gene, which questions would you ask and how would you proceed?
15. How would you convert a VCF file to the Plink binary format? How would you do the reverse, and what kind of problems do you anticipate?
16. Write a snippet to reformat a PED file so as to output a file with the following header `sample_name genotype_SNP1 genotype_SNP2 ...` where genotypes are coded VCF-style (e.g `A/C`, the order of the alleles in the output is not important).
17. A genetic association pipeline starts with a VCF and produces summary statistics for every position. The VCF contains multiallelics and indels. Unfortunately, a program in the pipeline trims all alleles to their first character. Why might allele frequencies not always be equal for a given variant? Find a way to correct the alleles in the association file by using the information from the VCF. Select columns are provided for [the association file](https://github.com/hmgu-itg/challenge/raw/master/data/association.txt.gz). We also provide [a file](https://github.com/hmgu-itg/challenge/raw/master/data/fromVCF.txt.gz) that was created from the VCF using `bcftools query -f '%CHROM %POS %REF %ALT %AN %AC\n'`.
18. [This file](https://github.com/hmgu-itg/challenge/raw/master/data/mock.eQTL.forChallenge.txt) contains eQTL overlap data for SNPs that arise as signals in GWAS for several phenotypes. Reformat this file to have one line per SNP/phenotype pair, and two additional columns formatted as such : `GENE1(tissue1, tissue2),GENE2(tissue1, tissue3)`, and `GENE1(2),GENE2(2)`. Each line should contain the SNP/phenotype pair, all genes found overlapping and their respective tissues, and all genes found overlapping with the number of tissues.
19. A researcher wants to conduct a disease association study. However, colleagues warn him that the dataset contains related individuals. He would like to remove relatedness in his dataset, but given his disease is rare, he would also like to maximise the number of cases kept in. Using [a list of samples with disease status](https://github.com/hmgu-itg/challenge/raw/master/data/relateds.pheno.tsv) and [a file containing pairs of individuals above a relatedness threshold](https://github.com/hmgu-itg/challenge/raw/master/data/relateds.tsv), create an exclusion list of samples to remove to help the researcher achieve their goal.

### Statistical genetics
1. You sample at random 10,000 variants from a deep (50x) whole-genome sequencing variant call file describing 1,000 individuals. What do you expect the distribution of minor allele frequency to look like? In particular, which minor allele counts are likely to be most frequent?
   <p>
    <ul>
    <b>Ans:</b> When you sample 10,000 variants at random from a deep (50×) whole‐genome sequencing dataset of 1,000 individuals and plot site frequency spectrum where Minor allele counts are on the X axis and Number of variants is on the Y, you would see more singletons (minor alleles found only once among all chromosomes). That's because singletons are expected to be the most frequent category because there are many recent or newly arisen mutations that haven’t had time to spread in the population. This essentially means the distribution should be highly skewed toward rare variants as site frequency spectrum 1/f law, which kind of means that the expected number of variants is proportional to 1/f. So, since singletons are rare i.e. has low frequencies they are pretty common and hence it is more likely to see them when you sample. Also, the human populations in particular have undergone recent expansions, amplifying the proportion of very rare variants.
    </ul>    
    </p>
3. You are running a single-point association study with a quantitative phenotype on the dataset above. Which filters, if any, would you apply to the dataset?
   <p>
    <p><b>Ans:</b> I will do the following QC steps:</p>
        <ul>
             <li>Call rate: In this step the goal is to clean the data by removing variants and individuals with too many missing genotypes.</li>
             <li>Depth or Coverage: Here, I will exclude genotypes that fall below a read-depth threshold (say 10x) and then remove variants with low call rates.</li>
             <li>Hardy-Weinberg equilibrium based filter: This step will perform HWE test in plink and remove spurious markers likely due to genotyping errors or non-random mating artifacts.</li>
             <li>Sex check: Here, I will use something like plink --check-sex to verify reported sex. The aim is to ensure there is no major mislabeled sexes in the data that could bias results, especially if you plan on including sex as a covariate or analyzing sex chromosomes.</li>
             <li>Minor allele frequency or MAF based filter: Possibly set a threshold like MAF ≥ 1% for a single‐marker association.</li>
             <li>Population structure / Relatedness: I will use PCA or relationship checks to catch outliers or duplicates.</li>
        </ul>
    </p>
    
5. A common practice when performing genetic association studies is to perform an ethnicity check as a quality control step. Can you explain how this is done?
    - You are dealing with whole-genome sequencing data in 2,326 Bulgarian samples. How would you perform such an ethnicity check, and which projection dataset would you use?
        	<ul>
	<b>Ans:</b> I will perform the following steps:
        <p>
    	<li> First, I will apply standard QC filters as described briefly in the second question's answer. This will allow me to extract a set of high-quality, common, and ideally independent SNPs from Bulgarian samples.</li> 
    
    	<li> Then I will look for Linkage Disequilibrium (LD) using plink in the data and try to prune. For this purpose, I'll use SNPs that are not in high LD with one another as this will help me get a clear picture of overall population structure.</li>
    
    	<li> Then I will try to merge my data with a reference dataset. For this pupose, first a reference dataset must be selected. The 1000 Genomes Project could be a good choice because it contains samples from multiple ethnic groups, including several European populations, e.g., CEU, GBR, FIN, IBS, TSI. For Bulgarian samples, I expect them to cluster with European reference populations. Then I'll combine the earlier pruned Bulgarian dataset with the reference panel. I need to make sure that to use the same set of SNPs and also the allele coding should be consistent across datasets.
    	</li>
    
    	<li>Then I will perform PCA on the merged dataset using plink. The analysis will generate principal components (PCs) that capture the major genetic variation.The PCA will yield coordinates for each individual and upon plotting these, there should be clusters corresponding to the known ethnic groups in the reference. I would expect the majority of the Bulgarian samples to cluster closely with European reference populations. I could also look for outliers or individuals (from Bulgarian dataset) fall far outside the European cluster, for example clustering with Africa or East Asian etc.These samples may need to be marked for further investigation or perhaps excluded from the study.
    	</li>
    	</ul>
         </p>
6. You perform a single-point association with blood lipids and find a variant with MAF=0.7% associated at p=1e-14. Do you expect the effect size to be large or small? What would be your next steps investigating this signal?
	<ul>
	<b>Ans:</b> When a variant with a minor allele frequency (MAF) of only 0.7% shows an extremely significant association (p=1e-14), this typically implies large effect size. That's because achieving p=1e-14 with so few minor alleles in the dataset usually requires a large per-allele effect, assuming the sample size is in the typical range for a genetic study (thousands to tens of thousands of individuals).Thsi also means, you likely have a relatively large cohort. Otherwise, you wouldn’t get such robust significance for a variant present in <1% of chromosomes.
	
	As a further investigation, I will try to do the following:

	<li>The most important next step would be to replicate the association in an independent sample set. If the signal reproduces, it strengthens confidence that this is a true finding rather than a statistical or technical artifact.</li>
	
	<li>I will also check the functional association of the findings, i.e., whether the variant lies in or near a known lipid-related gene or regulatory region. For this purspoe, I'll look up public databases such as Ensembl for annotation and predicted functional impact.</li>
	
	<li>I will suggest the experimentalists to perform orthogonal experiments (if feasible) to confirm a causal relationships.</li>
	
	<li>I’d look at the actual lipid measurements in carriers (LDL, HDL, etc.) to see if the difference is just statistically detectable or also clinically significant, for example, dangerously high LDL.</li>	
	</ul>
8. You are running an inverse-variance based single-point meta-analysis of the above dataset together with a UK study of 12,400 individuals. The variant above is found in the UK dataset, but its association is weak (1e-3). However, another variant located 1kb downstream is strongly associated with the same trait (p=1e-15) but weakly associated in your study (p=5e-4). Which of these 2 variants do you expect will have the strongest meta-analysis p-value? What do you think is happening in this region, how can you test it, and which model could you apply if it is the case?
9. An analyst studies a population of remote villages in Eastern Europe. They are interested in a particular variant, and compare the frequency in their villages (3.5%) to the EUR population frequency in the 1000 Genomes (0.03%). They conclude that the variant has increased in frequency in their villages. Do you agree, and if not, what would your advice be?
10.  The same analyst sends you association summary statistics for random glucose.
    - Which checks would you perform on such a dataset?
    - You wish to include this dataset in a meta-analysis. Which additional information should you ask for in your next email to your colleague?
    - In this dataset, you observe  &#955;=1.25. The analyst has adjusted for age, age-squared, sex, and fasting status. What would you suggest they do?
11. You are a co-author on a manuscript. You receive a draft, in which the main author has used the traditional &#945;=5e-8 as the significance threshold. The paper describes an analysis of 10 related blood phenotypes (platelet count, platelet volume, immature platelet fraction ...) using the fixed markers of the Infinium ImmunoArray on 897 individuals. What do you think about the chosen threshold, and what would you suggest to the first author? What would be your comments on the design of the study?
	<ul>
	<b>Ans:</b> I will let the first author know the following concerns:

	<li>For 10 phenotypes a stricter corerction (Bonferroni) would change the threshold to 5e-9. However, If the 10 phenotypes are correlated (e.g., platelet count and platelet volume might be strongly related), then we might not need a full 10-fold Bonferroni correction. Methods such as Sidak or matSpD that account for correlated phenotypes can yield an "effective number of independent tests" less than 10.</li>
	
	<li>Another point is that the Infinium ImmunoArray typically has a specific focus on immunologically relevant variants, and fewer total SNPs than a standard GWAS array which might have 500k–1M+ variants. Therefore, if we are analyzing fewer variants, for example, ~200k or 300k, then the multiple testing burden should be lesser than a million. Therefore, I might argue about using a slightly less stringent threshold than 5e-8, especially if we only test the variants actually on that array.</li>
	
	<li>The other concern that I have is about the small sample size. With fewer than 1,000 individuals, the statistical power to detect common variant associations, especially at a threshold of 5e-8 will be very limited. I think, for rare variants or small effects, it is extremely challenging to reach 5e-8 significance in only 897 subjects. Therefore, it may be possible that reviewers of the manuscript ask for replication or meta-analysis with larger samples.</li>
