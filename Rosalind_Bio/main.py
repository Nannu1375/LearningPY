from customTool import *
def main():

    #read the file to get sequence
    file_dir = "/media/nannu/Briefcase/Python_wrkspc/rosalind/Seq_Analysis/rosalind_gc.txt"

    result_id, highgccount = find_highest_gc_content(file_dir)

    print(f"{result_id}\n{highgccount:.6f}%")



    # with open(file_dir, 'r') as file:
    #     seq = file.readline()

    # #count the bases and print them
    # rComp = reverse_comp(seq)
    # print(rComp)

    # with open(file_dir, 'r') as file:
        # months, pairs = map(int, file.readline().strip().split())
    # 
    # print(modFibonacci(months, pairs))





#call the main function
if __name__ == "__main__":
    main()
