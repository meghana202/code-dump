#Converts a number into its nitrogenous base (DNA)
function NumToChar(Num)
    Char = Dict(0=>'A', 1=>'C', 2=>'G',3=>'T')
    return Char[Num]
end
NumToChar(1)
#Output: 'C': ASCII/Unicode U+0043 (category Lu: Letter, uppercase)

#Converts a character (DNA base) into a number
function CharToNum(Char)
    BaseFour = Dict('A'=>0,'C'=>1,'G'=>2,'T'=>3)
    return BaseFour[Char]
end
CharToNum('C')
#Output: 1

#Converts a string of **DNA** to base 4
function PatternToNumber(Pattern)
    if length(Pattern) == 0
        return 0
    end
    symbol = Pattern[end]
    Prefix = Pattern[1:end-1]
    answer = 4*PatternToNumber(Prefix) + CharToNum(symbol)
    return answer
end
PatternToNumber("CATGAGGAAACCGGACTTAA")
#Output: 335678906864

#Converts a number(index) of length k back into DNA bases
function NumberToPattern(index,k)
    if k == 1
        return NumToChar(index)
    end
    prefixIndex = floor(index/4)
    r = index%4
    symbol = NumToChar(r)
    PrefixPattern = NumberToPattern(prefixIndex, k-1)
    return PrefixPattern * symbol
end
NumberToPattern(335678906864,20)
#Output: "CATGAGGAAACCGGACTTAA"

#Finds the frequency of a kmer in the text (in numerical order eg. FrequencyArray[1] corresponds to AA,[2] is AC and so on)
function ComputingFrequencies(Text,k)
    FrequencyArray = zeros(Int64, 4^k)
    for i = 1:(length(Text)-(k-1))
        Pattern = Text[i:i+k-1]
        j = PatternToNumber(Pattern)
        FrequencyArray[j+1] += 1
    end
    return FrequencyArray
end
ComputingFrequencies("ACGCGGCTCTGAAA",2)
#Output: 16-element Vector{Int64}:
 2
 1
 0
 0
 0
 0
 2
 2
 1
 2
 1
 0
 0
 1
 1
 0


#location of least skew in the genome of a bacteria
function Skew(genome)
    skew_i = 0
    skew = []
    for i = 1:length(genome)
        if genome[i]=='C'
            skew_i -= 1
        elseif genome[i]=='G'
            skew_i+=1
        end
        skew = append!(skew,skew_i)
    end
    x = minimum(skew)
    for i = 1:length(skew)
        if skew[i]==x
            return i
        end
    end
end
#genome = open("Salmonella_enterica.txt", "r")
genome = open("ecoli.txt","r")
Genome = read(genome,String)
newline = "\n"
Genome = replace(Genome,newline => "" )
Skew(Genome)
#Genome[3923620:3923620+500] <-- this is the window for OriC
#Output: 3923620

#The number of different base pairs in sequences of the same length (p and q)
function HammingDistance(p,q)
    HammingDistance = 0
    for i = 1:length(p)
        if p[i] != q[i]
            HammingDistance += 1
        end
    end
    return HammingDistance
end
HammingDistance("AGCTAGAAAATCGAT","CCCTAGAAAATCGAT")
#Output: 2

#returns the reverse complement of DNA string
function ReverseComplement(Pattern)
    complement = []
    x = length(Pattern)
    for i = 0:x-1
        if string(Pattern[x-i])=="A"
            append!(complement, "T")
        elseif string(Pattern[x-i])=="T"
            append!(complement,"A")
        elseif string(Pattern[x-i])=="C"
            append!(complement,"G")
        elseif string(Pattern[x-i])=="G"
            append!(complement, "C")
        end
    end
return join(complement)
end
ReverseComplement("TTATCCACA")
#Output: "TGTGGATAA"

#a list of all 'neighbor' strings with hamming distance <= d
function Neighbors(Pattern,d)
    if d == 0
        return Pattern
    end
    if length(Pattern) == 1
        return "A C G T"
    end
    Neighborhood = String[]
    pattern = Pattern[2:length(Pattern)]
    SuffixNeighbors = Neighbors(pattern,d)
    nucleotide = String["A","C","G","T"]
    for Text in SuffixNeighbors
        if HammingDistance(pattern,string(Text))<d
            for i = 1:4
                text = nucleotide[i]*Text
                if string(text[length(text)])!=" "
                    push!(Neighborhood,string(text))
                end
            end
        else
            Text = Pattern[1]*Text
            if string(Text[length(Text)])!= " "
                push!(Neighborhood,string(Text))
            end
        end
    end
    return Neighborhood
end 
Neighbors("CCAGTCAATG",1) 
#Output: 31-element Vector{String}:
 "CCAGTCAATA"
 "CCAGTCAATC"
 "CCAGTCAAAG"
 "CCAGTCAACG"
 "CCAGTCAAGG"
 "CCAGTAAATG"
 "CCAGACAATG"
 "CCAGCCAATG"
 "CCAGGCAATG"
 "CCAATCAATG"
 "CCACTCAATG"
 "CAAGTCAATG"
 "ACAGTCAATG"
 ⋮
 "CCGGTCAATG"
 "CCTGTCAATG"
 "CCATTCAATG"
 "CCAGTGAATG"
 "CCAGTTAATG"
 "CCAGTCCATG"
 "CCAGTCGATG"
 "CCAGTCTATG"
 "CCAGTCACTG"
 "CCAGTCAGTG"
 "CCAGTCATTG"
 "CCAGTCAATT"

#kmers with most neighbors in text, kmer doesn't necessarily appear in text, has hamming distance <= d
function FrequentMismatch(Text, k, d)
    FrequencyArray = zeros(Int64,4^k)
    for i = 1:(length(Text)-(k-1))
        Pattern = Text[i:i+k-1]
        Neighborhood = Neighbors(Pattern,d)
        for ApproximatePattern in Neighborhood
            j = PatternToNumber(ApproximatePattern)
            h = PatternToNumber(ReverseComplement(ApproximatePattern))
            FrequencyArray[j+1] += 1
            FrequencyArray[h+1] += 1
        end
    end
    maxCount = maximum(FrequencyArray)
    Mismatches = String[]
    for i = 1:length(FrequencyArray)
        if FrequencyArray[i]==maxCount
            push!(Mismatches, NumberToPattern(i,k))
        end
    end
    return Mismatches
end
FrequentMismatch("TGCAATGTCAAATACATACCAGAAAAGAAAGGCAAGTAAAGAGTATGGTTTATATAATCAATGTAAGAAACTAAATGATGATGAATTATTTCGCTTACTTGATGATCACAATTCCTTGAAAAGGATTTCATCTGCCAGAGTATTACAGTTAAGAGGTGGGCAAGACGCTGTTAGATTGGCAATTGAGTTCTGCTCTGATAAAAATTATATCCGTAGAGATATCGGAGCATTTATACTCGGGCAAATAAAAATTTGCAAAAAATGCGAAGATAATGTTTTTAATATTTTGAACAATATGGCATTGAATGATAAGAGCGCTTGCGTTCGAGCTACGGCAATCGAGTCAACGGCTCAGCGATGCAAGAAAAACCCAATTTATTCACCTAAAATAGTAGAACAATCTCAAATTACTGCTTTTGATAAATCGACTAATGTCAGACGTGCTACAGCATTTGCTATTTCTGTTATCAATGATAAAGCAACAATTCCACTATTGATTAA",9,1)
#Output: 10-element Vector{String}:
 "AATAATTTT"
 "AATATCAAC"
 "ATAGCAGAC"
 "ATTGATAAC"
 "CAAATTCAA"
 "TGATCAAAC"
 "TTCTGCTCA"
 "TTTATCACA"
 "TTTGATCAA"
 "TTTTGATCC"
          
#checks what kmers from FrequentMismatch(Text) actually appear in Text 
function Check(Text)
    Array = String[]
    Pattern = FrequentMismatch(Text, 9, 1)
    print(Pattern)
    for i = 1:length(Pattern)
        if occursin(Pattern[i], Text)
            push!(Array, Pattern[i])
        end
    end
    return Array
end
Check("CAATGATGATGACGTCAAAAGGATCCGGATAAAACATGGTGATTGCCTCGCATAACGCGGTATGAAAATGGATTGAAGCCCGGGCCGTGGATTCTACTCAACTTTGTCGGCTTGAGAAAGACCTGGGATCCTGGGTATTAAAAAGAAGATCTATTTATTTAGAGATCTGTTCTATTGTGATCTCTTATTAGGATCGCACTGCCCTGTGGATAACAAGGATCCGGCTTTTAAGATCAACAACCTGGAAAGGATCATTAACTGTGAATGATCGGTGATCCTGGACCGTATAAGCTGGGATCAGAATGAGGGGTTATACACAACTCAAAAACTGAACAACAGTTGTTCTTTGGATAACTACCGGTTGATCCAAGCTTCCTGACAGAGTTATCCACAGTAGATCGCACGATCTGTATACTTATTTGAGTAAATTAACCCACGATCCCAGCCATTCTTCTGCCGGATCTTCCGGAATGTCGTGATCAAGAATGTTGATCTTCAGTG")
#Output: 1-element Vector{String}:
 "CTGGGATCC"    


