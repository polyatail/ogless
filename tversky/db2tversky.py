# -*- coding: utf-8 -*-


#    tversky_matrix = numpy.zeros((len(index["genomes"]), len(index["genomes"])))
#
#    for nameA, nameB in itertools.combinations(index["genomes"], 2):
#        nameA_num = index["genome_to_num"][nameA]
#        nameB_num = index["genome_to_num"][nameB]
#        
#        sizeA = len(index["genome_to_cds"][nameA_num])
#        sizeB = len(index["genome_to_cds"][nameB_num])
#
#        hitsA = index["blast_hits"][nameA_num][nameB_num]
#        hitsB = index["blast_hits"][nameB_num][nameA_num]
#
#        if hitsA > sizeA or \
#           hitsB > sizeB:
#            raise ValueError("shared genes cannot exceed genome size")
#
#        shared_genes = min([hitsA, hitsB])
#
#        a_not_b = sizeA - shared_genes
#        b_not_a = sizeB - shared_genes
#        
#        tversky50 = float(shared_genes) / float(shared_genes + 0.5 * a_not_b + 0.5 * b_not_a)
#
#        tversky_matrix[nameA_num][nameB_num] = (1 - tversky50) ** 3
#        tversky_matrix[nameB_num][nameA_num] = (1 - tversky50) ** 3
#
#        #if tversky50 < 0.4: continue
#        #print "\t".join(map(str, [nameA, "(pp)", nameB, tversky50]))
#
#    print "\t%s" % len(index["genomes"])
#
#    for genome_num, genome_name in enumerate(index["genomes"]):
#        row = [genome_name] + [x[:8] for x in map(str, tversky_matrix[genome_num])]
#        print "\t".join(row)