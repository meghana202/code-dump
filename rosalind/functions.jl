# housekeeping
using ArgParse
using CodecZlib
using FASTX
using BioSequences
using HTTP

cli = ArgParseSettings()
@add_arg_table cli begin
        "--file"
                help = "file"
        "--fasta"
                help = "fasta"
end
arg = parse_args(ARGS, cli)

#= READER
reader = nothing
if splitext(arg["fasta"])[2] == ".gz"
	reader = FASTA.Reader(GzipDecompressorStream(open(arg["fasta"])))
else
	reader = FASTA.Reader(open(arg["f"]))
end
=#
#my own fasta thing
function fasta_decoder(file)
	file = readlines(open(file, "r"))
	id = []
	seqs = []
	i = 1
	while i <= length(file)
		seq = ""
		if file[i][1] == '>' 
			push!(id, file[i])
		end
		i += 1
		while file[i][1] != '>'
			seq = seq*file[i]
			i += 1
			if i > length(file) break end
		end
		push!(seqs, seq)
	end
	return(id, seqs)
end
#println(fasta_decoder("test.fa"))
#end, rosalind functions below


#INI5, file
function second_line()
    file = open(arg["file"], "r")
    lines = readlines(file)
    #print(typeof(lines))
    for i = 1:length(lines)
    	if i%2 == 0 
    		println(lines[i])
    	end
    end
end
#second_line()

#INI6, file
function word_dict()
	str = open(arg["file"], "r")
	Str = read(str, String)
	array_ = split(Str, ('\n', ' '))
	answer = Dict()
	for word in array_
		if haskey(answer, word) 
			answer[word] += 1
		elseif word == ""
		#do nothing
		else answer[word] = 1
		end
	end	
	for (key, value) in answer
		println("$key $value")
	end
end
#word_dict()

#DNA, string
function dna_count()
	str = open(arg["file"], "r")
	Str = read(str, String)
	answer = Dict()
	for nt in Str
		if haskey(answer, nt)
			answer[nt]+=1
		elseif nt == '\n'
		else answer[nt] =1
		end 
	end
	print("$(answer['A']) $(answer['C']) $(answer['G']) $(answer['T']) ")
end
#dna_count()

#RNA, file
function torna()
	str = open(arg["file"], "r")
	Str = read(str, String)
	answer = ""
	for letter in Str
		if letter != 'T'
			answer*=letter
		elseif letter == 'T'
			answer*='U'
		end
	end
	println(answer)
end
#torna()

#REVC, file
function complement()
	str = open(arg["file"], "r")
	Str = read(str, String)
	answer = ""
	for i = 0:length(Str)-1
		if Str[length(Str)-i] == 'T'
			answer*='A'
		elseif Str[length(Str)-i] == 'A'
			answer*='T'
		elseif Str[length(Str)-i] =='C'
			answer*='G'
		elseif Str[length(Str)-i] =='G'
			answer*='C'
		end
	end
	println(answer)
end
#complement()

#FIB, nothing
function wabbits(n, k)
	if n== 1
		wabbit = 1
		return wabbit
	elseif n == 2
		wabbit = 1
		return wabbit
	elseif n != 1
		for i = 2:n
			wabbit = wabbits(n-1,k)+ wabbits(n-2,k)*k
			return (wabbit)
		end
	end
end
#println(wabbits(35,4))

#GC, fasta, gc_ finds the content, main returns the answer from a fasta file, reader
function gc_(seq)
	content = 0
	for i = 1:length(seq)
		if seq[i] == 'C'||seq[i] == 'G' content += 1 end
	end
	answer = (content/length(seq))*100
	return answer
end

function main()
	gc = []
	id = []
	for record in reader
		rec = FASTA.identifier(record)
		dna_seq = FASTA.sequence(String, record)
		push!(gc, gc_(dna_seq))
		push!(id, rec)
	end
	println("$(id[argmax(gc)])")
	println("$(maximum(gc))")
end
#main()

#HAMM, file
function hamm()
	file = open(arg["file"], "r")
	File = read(file, String)
	strings = split(File, '\n')
	string1 = strings[1]
	string2 = strings[2]
	distance = 0
	for i = 1:length(string1)
		if string1[i] != string2[i] distance += 1 end
	end
	println(distance)
end
#hamm()

#IPRB, nothing
function mendel(k, m, n)
	total = (k+n+m)*(k+n+m-1)
	k_sum = k*(k-1) + k*m + k*n 
	n_sum = 0*n*(n-1) + n*k + .5*n*m
	m_sum = .75*m*(m-1) + m*k + .5*m*n
	kmn_sum = k_sum + n_sum + m_sum
	answer= kmn_sum/total
	println(answer)
end
#mendel(19, 20, 16)

#PROT, dict amino (see below) needed 
#= 
amino = Dict{Any, Any}("GAC" => "D", 
"UAU" => "Y", "UUA" => "L", "CUU" => "L", "CCU" => "P", 
"GGC" => "G", "GCU" => "A", "CGG" => "R", "CAG" => "Q", 
"ACG" => "T", "GAU" => "D", "UGU" => "C", "AGA" => "R", 
"GUU" => "V", "CUA" => "L", "UUG" => "L", "CUC" => "L", 
"CGA" => "R", "CAA" => "Q", "UUU" => "F", "UUC" => "F", 
"AUU" => "I", "CUG" => "L", "UAA" => "Stop", "AAC" => "N", 
"GUC" => "V", "UCC" => "S", "CCA" => "P", "CAC" => "H", 
"GCA" => "A", "CCG" => "P", "GAA" => "E", "GUA" => "V", 
"ACA" => "T", "GUG" => "V", "UCU" => "S", "AGG" => "R", 
"AUC" => "I", "AUG" => "M", "AUA" => "I", "UCG" => "S", 
"AGC" => "S", "ACU" => "T", "GCC" => "A", "AAG" => "K", 
"GGU" => "G", "AAU" => "N", "AAA" => "K", "UAG" => "Stop", 
"UGC" => "C", "ACC" => "T", "UGA" => "Stop", "AGU" => "S", 
"CCC" => "P", "UAC" => "Y", "GAG" => "E", "CGU" => "R", 
"GCG" => "A", "CAU" => "H", "CGC" => "R", "UGG" => "W",
"GGA" => "G", "UCA" => "S", "GGG" => "G")
=#
function protein()
	file = open(arg["fasta"], "r")
	File = read(file, String)
	answer = ""
	for i = 1:3:length(File)-2
		codon = File[i:i+2]
		if amino[codon] != "Stop" answer*=(amino[codon]) end
	end
	println(answer)
end
#protein()

#SUBS, fasta
function motifs()
	file = open(arg["fasta"], "r")
	File = read(file, String)
	File = split(File, '\n')
	seq = File[1]
	motif = File[2]
	answer = []
	for i = 1:length(seq)-length(motif)+1
		#println(seq[i:i+length(motif)-1])
		if seq[i:i+length(motif)-1] == motif push!(answer, i)end
	end
	for number in answer print("$number ") end
end
#motifs()

#CONS, uses array file (below), reader
#=
file = []
for record in reader
	rec = FASTA.identifier(record)
	dna_seq = FASTA.sequence(String, record)
	push!(file, rec)
	push!(file, dna_seq)
end 
=#
function consensus()
	profile_ = zeros(Int8, 4, length(file[2]))
	for i = 2:2:length(file)
		for j = 1:length(file[i])
			str = file[i]
			if str[j] == 'A'
				profile_[1,j] += 1
			elseif str[j] == 'C'
				profile_[2,j] += 1
			elseif str[j] == 'G'
				profile_[3,j] += 1
			elseif str[j] == 'T'
				profile_[4,j] += 1
			end
		end
	end
	for col in eachcol(profile_)
		loc = (argmax(col))
		if loc == 1
			print('A')
		elseif loc == 2
			print('C')
		elseif loc == 3
			print('G')
		elseif loc == 4
			print('T')
		end
	end
	println()
	x = ["A: ", "C: ", "G: ", "T: "]
	for i = 1:4
		print(x[i])
		for j = 1:length(profile_[i, :])
			print("$(profile_[i,j]) ")
		end
		println()
	end
end
#consensus()

#GRPH, uses Dict file, see below, reader
#= 
file = Dict()
for record in reader
	rec = FASTA.identifier(record)
	dna_seq = FASTA.sequence(String, record)
	file[rec] = (dna_seq[1:3], dna_seq[end-2:end])
end
=#
function graphs()
	for (key,value) in file
		first_ = value[2]
		for (key1,value1) in file
			last_ = value1[1]
			if first_ == last_ && key != key1
				#println(first_*" "*last_)
				println("$(key) $(key1)")
			end
		end
	end
end
#graphs()

#IEV, nothing
function E(x)
	answer = 2*(sum(x[1:3]) + x[4]*.75	 + x[5]*.5)
	println(answer)
end
#E([19918 19099 19809 18467 18778 16436])

#LCSM, uses array seq(see below), reader
#=
seq = []
for record in reader
	rec = FASTA.identifier(record)
	dna_seq = FASTA.sequence(String, record)
	push!(seq, dna_seq)
end
=#

function shared_motif()
	answers = []
	sequence = seq[1]
	for i = 0:length(sequence)-1
		x = length(sequence)-i
		for j = 1:length(sequence)-x
			seq_ = sequence[j:j+x]
			#println(seq_)
			count_ = 0
			for i = 1:length(seq) if occursin(seq_, seq[i]) count_+=1  end end
			if count_ == length(seq) push!(answers, (length(seq_), seq_)) end
		end
	end 
	println(sort(answers)[(length(answers))][2])
end
#shared_motif()

#ORF, fasta, reader inside function
function orf()
	reader = nothing
	if splitext(arg["fasta"])[2] == ".gz"
		reader = FASTA.Reader(GzipDecompressorStream(open(arg["fasta"])))
	else
		reader = FASTA.Reader(open(arg["fasta"]))
	end
	orfs = []
	for record in reader
		#println(">", FASTA.identifier(record))
		dna_seq = FASTA.sequence(String, record)
		comp = complement(dna_seq)
		for i = 1:length(dna_seq)-2
			codon = dna_seq[i:i+2]
			if codon == "ATG"
				for j = i:3:length(dna_seq)-2
					newcod = dna_seq[j:j+2]
					if newcod == "TAA" || newcod == "TAG" || newcod == "TGA" 
						orf = dna_seq[i:j-1]
						push!(orfs, orf)
						break
					end 
				end 
			end 
			codon2 = comp[i:i+2]
			if codon2 == "ATG"
				for j = i:3:length(comp)-2
					newcod = comp[j:j+2]
					if newcod == "TAA" || newcod == "TAG" || newcod == "TGA" 
						orf = comp[i:j-1]
						push!(orfs, orf)
						break
					end 
				end 
			end 
		end
	end
	orfs = Set(orfs)
	for orf in orfs
		orf2 = ""
		for letter in orf
			if letter == 'T'
				orf2*='U'
			else
				orf2*= letter
			end
		end
		for i = 1:3:length(orf)
			print("$(amino[orf2[i:i+2]])")
		end
		println()
	end
end
#orf()

#PRTM, file, uses weight_dict below
#=
weight_dict = Dict{Any, Any}("Q" => 128.05858, "W" => 186.07931, 
"T" => 101.04768, "C" => 103.00919, "P" => 97.05276, 
"V" => 99.06841, "L" => 113.08406, "M" => 131.04049, 
"N" => 114.04293, "H" => 137.05891, "A" => 71.03711, 
"D" => 115.02694, "G" => 57.02146, "E" => 129.04259, 
"Y" => 163.06333, "I" => 113.08406, "S" => 87.03203, 
"K" => 128.09496, "R" => 156.10111, "F" => 147.06841)
=#
function weight()
	file = open(arg["file"], "r")
	File = read(file, String)
	File = split(File, ('\n'))
	File = join(File)
	sum = 0 
	for i = 1:length(File)
		sum += weight_dict[string(File[i])]
	end
	println(sum)
end
#weight()

#MPRT, file
function protein_motif()
	file = read(open(arg["file"], "r"), String)
	file = split(file, "\n")
	#used length(file)-1 because the dataset comes with an empty line
	file = file[1:length(file)-1]
	website = "https://www.uniprot.org/uniprot/"
	for ID in file 
		r = HTTP.open("GET", website*ID*".fasta") do io
   			while !eof(io)
   				x = String(readavailable(io))
   				x = split(x, "\n")
   				x = join(x[2:length(x)-1])
   				locations = [ID]
   				for i = 1:length(x)-3
   					protein_seq = x[i:i+3]
   					if protein_seq[1] == 'N' && protein_seq[2] != 'P' && (protein_seq[3] == 'S'||protein_seq[3]=='T') &&protein_seq[4]!=4
   						push!(locations, string(i))
   					end
   				end 
   				if length(locations) > 1
   					println(locations[1])
   					for i = 2:length(locations) print("$(locations[i]) ") end 
   					println()
   				end
   			end
		end
	end
end
#protein_motif()

#MRNA, dict aa(see below), file
#=
Dict{Any, Any}("Q" => 2, "W" => 1, "T" => 4, 
"P" => 4, "C" => 2, "V" => 4, "L" => 6, 
"M" => 1, "Stop" => 3, "N" => 2, "H" => 2, 
"A" => 4, "D" => 2, "G" => 4, "E" => 2, "Y" => 2,
"I" => 3, "S" => 6, "K" => 2, "R" => 6, "F" => 2)
=#
function mRNA()
	count_ = 3
	file = read(open(arg["file"], "r"), String)
	file = join(split(file, '\n'))
	begin_ = findfirst(isequal('M'), file)
	end_ = findfirst(isequal("Stop"), file)
	if end_ != nothing
		protein = file[begin_:end_-1]
	elseif end_ == nothing
		protein = file[begin_:end]
	end
	for letter in protein
		count_ *= aa[string(letter)]
		mod_ = count_%1000000
		count_ = mod_
	end
	println(count_)
end
#mRNA()

#REVP (rc and REVP), uses file
function rc(dna)
	rc_ = ""
	for i = 0:length(dna)-1
		if dna[length(dna)-i] == 'A' rc_ *= 'T'
		elseif dna[length(dna)-i] == 'T' rc_ *= 'A'
		elseif dna[length(dna)-i] == 'G' rc_ *= 'C'
		elseif dna[length(dna)-i] == 'C' rc_ *= 'G' end
	end
	return rc_
end 

function REVP()
	file = read(open(arg["file"], "r"), String)
	file = join(split(file,'\n')[2:end])
	for i = 1:length(file)-11
		for j = 4:2:12
			if file[i:i+j-1] == rc(file[i:i+j-1]) println("$i $j") end
		end
	end
	last_ = file[length(file)-10:end]
	for i = 1:11
		for j = 4:2:12-i
			if last_[i:i+j-1] == rc(last_[i:i+j-1]) println("$(length(file)-(11-i)) $j") end
		end
	end
		
end
#REVP()

#SPLC, uses file, fasta_decoder (above), dict amino (see #PROT)
function splicing()
	file = fasta_decoder(arg["file"])
	dna = file[2][1]
	introns = file[2][2:end]
	rna = ""
	for intron in introns
		dna = replace(dna, intron => "")
	end
	for letter in dna
		if letter == 'T' rna = rna*'U'
		else rna = rna*letter
		end
	end
	for i = 1:3:length(rna)-2
		if amino[rna[i:i+2]] == "Stop" break end
		print(amino[rna[i:i+2]])
	end
	println()
end
#splicing()
		
#PERM 
function perm(n,m)
	list = []
	nums = 1:n
	if m == 1
		for i = 1:n push!(list, string(i)) end
		return list
	elseif m > 1
		temp = perm(n,m-1)
		for i = 1:length(temp)
			for j in nums
				if !occursin(string(j), temp[i])
					push!(list, temp[i]*" "*string(j)) 
				end 
			end
		end 
		return list
	end			
end


function output_formatting(n)
	println(factorial(n))
	x = perm(n,n)
	println(length(x))
	for i = 1:factorial(n) println(x[i]) end
end
output_formatting(7)


			
 






		