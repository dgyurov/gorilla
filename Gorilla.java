import java.io.File;
import java.io.IOException;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Scanner;

public class Gorilla {
	
	private static HashMap<String, Integer> costMap;
	private static int[][] M;

	public static void main(String[] args){
		//Set the data folder
		Path dir = Paths.get(".");
        
		//Initialize a list of strings
		List<String> files = new ArrayList<String>();
        
		//Get all file paths inside the directory and populate the list
        files = getAllFileNames(files, dir);

		//Generate a cost map from the BLOSUM matrix
		costMap = CostMatrix.generateCostMap();
        
		//Run the algorithm for all files
        for (String file : files) {
        	algorithm(file);
		}
	}
	
	public static List<Organism> parsefile(String fileName) {
		
		//Create a new File object using the passed filename
		File file = new File(fileName);
		
		//Initialize an arraylist of organisms
		ArrayList<Organism> organisms = new ArrayList<Organism>();
		Scanner scanner = null;

		//Initialize the scanner
		try {
			scanner = new Scanner(file);
		} catch (Exception e) {
			System.exit(0);
		}
		try {
			//Get the first line in the file
			String line = scanner.nextLine().trim();
			
			while (scanner.hasNext()) {
				
				//Split by any kind of white space
				String[] splitLine = line.split("\\s+");
				
				//Get the name from the line
				String name = splitLine[0].substring(1);
				
				//Initialize a string which will hold the protein sequence
				String proteinSequence = "";
				while (scanner.hasNext()) {
					line = scanner.nextLine().trim();
					
					//If the next line is a new organism
					if (line.startsWith(">")) {
						break;
						
					//Otherwise add the line to the proteinSequence for this organism
					} else {
						proteinSequence += line;
					}
				}
				
				//Create a new organism object and add it to the list of organisms
				organisms.add(new Organism(name, proteinSequence.toCharArray()));
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		return organisms;

	}
	
	
	
    public static List<String> getAllFileNames(List<String> fileNames, Path dir) {
        try {
        	
            //Create a new dir stream
        	DirectoryStream<Path> stream = Files.newDirectoryStream(dir);
            
            for (Path path : stream) {
            	//If its a sub dir, call the method recursively 
                if (path.toFile().isDirectory()) {
                    getAllFileNames(fileNames, path);
                } else {
                	//Otherwise add the file path to the list if it's a FASTA input one
                	String fileName = path.toAbsolutePath().toString();
                	if(fileName.contains("FASTAs-in")){
                        fileNames.add(fileName);
                	}
                }
            }
            stream.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return fileNames;
    }
	
	public static void algorithm(String file){
		//Parse the current file and populate the List of organisms
		List<Organism> organisms = parsefile(file);

		//Align an organism with the rest of the organisms in the file
		for (int i = 0; i < organisms.size() - 1; i++) {
			Organism organism1 = organisms.get(i);
			int length1 = organism1.sequence.length;
			
			for (int j = i + 1; j < organisms.size(); j++) {
				Organism organism2 = organisms.get(j);
				int length2 = organism2.sequence.length;

				//Initialize the M, which will hold the optimal score
				M = new int[length1][length2];
				for (int k = 0; k < length1; k++) {
					for (int l = 0; l < length2; l++) {
						M[k][l] = Integer.MIN_VALUE;
					}
				}
				
				//Calculate the most optimal score based on the BLOSUM matrix
				int opt = opt(length1 - 1, length2 - 1, organism1.sequence, organism2.sequence);
				
				//Print the result
				System.out.println(organism1.name + "--" + organism2.name + ": " + opt);
				
				//Print the optimal sequences
				System.out.println(printResult(organism1.sequence,organism2.sequence));
			}

		}
		
	}

	//Based on Needleman–Wunsch algorithm
	public static String printResult(char[] sequence1, char[] sequence2) {

		//Initialize the results
		String organism1 = "", organism2 = "";
	
		//Initialize the pointers of the sequences, start from the last protein
		int s1Pointer = sequence1.length-1, s2Pointer = sequence2.length-1;	
		
		//While there are more than 1 proteins left in either sequence
		while(s1Pointer>0 || s2Pointer>0){
			
			//If the best alignment is using the proteins from both sequences
			if (s1Pointer > 0 && s2Pointer > 0 && M[s1Pointer][s2Pointer] == M[s1Pointer - 1][s2Pointer - 1] + costMap.get(String.valueOf(sequence1[s1Pointer]) + String.valueOf(sequence2[s2Pointer]))) {
				organism1 = String.valueOf(sequence1[s1Pointer--]) + organism1;
				organism2 = String.valueOf(sequence2[s2Pointer--]) + organism2;

			//If the best alignment is using a dash for the second sequence
			} else if (s1Pointer > 0 && M[s1Pointer][s2Pointer] == M[s1Pointer - 1][s2Pointer] + costMap.get("A*")) {
				organism1 = String.valueOf(sequence1[s1Pointer--]) + organism1;
				organism2 = "-" + organism2;
			
			//If the best alignment is using a dash for the first sequence	
			} else if (s2Pointer > 0 && M[s1Pointer][s2Pointer] == M[s1Pointer][s2Pointer - 1] + costMap.get("A*")) {
				organism1 = "-" + organism1;
				organism2 = String.valueOf(sequence2[s2Pointer--]) + organism2;
			
			//Otherwise break
			} else {
				break;
			}
		}
		
		//add the first protein to the resulting sequence
		organism1 = sequence1[0] + organism1;
		organism2 = sequence2[0] + organism2;
		
		return organism1 + "\n" + organism2;
	}

	//The bigger the score, the more close the organisms
	public static int opt(int n, int m, char[] seq1, char[] seq2) {
		
		//Base case if n is smaller than 0
		if (n < 0) {
			return (m + 1) * costMap.get("A*").intValue();
		
		//Base case if m is smaller than 0
		} else if (m < 0) {
			return (n + 1) * costMap.get("A*").intValue();
		}

		//If the M(n,m) was already computed, return it
		if (M[n][m] != Integer.MIN_VALUE) {
			return M[n][m];
		}

		// (i) (m,n) -> M
		int one = costMap.get(String.valueOf(seq1[n]) + String.valueOf(seq2[m])) + opt(n - 1, m - 1, seq1, seq2);
		
		// (ii) The n'th position is not matched
		int two = costMap.get(String.valueOf(seq1[n]) + "*") + opt(n - 1, m, seq1, seq2);
		
		// (iii) The m'th position is not matched
		int three = costMap.get(seq2[m]+"*") + opt(n, m - 1, seq1, seq2);
		
		//Get the maximum score out from the three
		int max = getMax(one, two, three);
		
		//Fill in the optimal score in M and return it
		M[n][m] = max;
		return max;
	}
	
	//Returns the maximum of three numbers
	public static int getMax(int a, int b, int c){
		int biggest = a;
		if (biggest < b) biggest = b;
		if (biggest < c) biggest = c;
		return biggest;
	}
}

class Organism {
	
	//An organism object holds the name and the sequence of proteins
	public String name;
	public char[] sequence;
	
	public Organism(String name, char[] sequence){
		this.name = name;
		this.sequence = sequence;
	}
}

/**
 * #  Matrix made by matblas from blosum62.iij
 * #  * column uses minimum score
 * #  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
 * #  Blocks Database = /data/blocks_5.0/blocks.dat
 * #  Cluster Percentage: >= 62
 * #  Entropy =   0.6979, Expected =  -0.5209
 */

class CostMatrix {

	//Populates a HashMap with all possible combinations, for faster request calls later on
	public static HashMap<String, Integer> generateCostMap(){
		HashMap<String, Integer> costMap = new HashMap<String, Integer>();
		
			for (int i = 0; i < matrix.length; i++) {
				for (int j = 0; j < matrix[0].length; j++) {
					String combination = String.valueOf(legend[i]) + String.valueOf(legend[j]);
					costMap.put(combination, matrix[i][j]);
				}
			}
		
		return costMap;
	}
	
	//All possible proteins in the matrix
	private static char[] legend = new char[] {
			'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*'
	};
	
	//The matrix represented as an array of arrays
	private static int[][] matrix = new int[][] {
		
		 /*     A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   * */
		/*A*/ { 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4 },
		/*R*/ {-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4 },
		/*N*/ {-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4 },
		/*D*/ {-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4 },
		/*C*/ { 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4 },
		/*Q*/ {-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4 },
		/*E*/ {-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4 },
		/*G*/ { 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4 },
		/*H*/ {-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4 },
		/*I*/ {-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4 },
		/*L*/ {-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4 },
		/*K*/ {-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4 },
		/*M*/ {-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4 },
		/*F*/ {-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4 },
		/*P*/ {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4 },
		/*S*/ { 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4 },
		/*T*/ { 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4 },
		/*W*/ {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4 },
		/*Y*/ {-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4 },
		/*V*/ { 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4 },
		/*B*/ {-2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4 },
		/*Z*/ {-1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4 },
		/*X*/ { 0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4 },
	  /* * */ {-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1 },	
		};
	}
