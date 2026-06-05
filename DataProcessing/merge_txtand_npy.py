import numpy as np
import argparse


def add_txt_to_npy(txt_folder, npy_folder, output_folder):
    # Define file paths
    #generate positive samples with location
    txt_files = [f"positive_chr{chrom}.txt" for chrom in range(1, 23)] + ["positive_chrX.txt"]
    npy_files = [f"KR_5kb_matrix_chr{chrom}_positive.npy" for chrom in range(1, 23)] + [
    # generate negative samples with location
    # txt_files = [f"negative_chr{chrom}-5KB-sort.txt" for chrom in range(1, 23)] + ["negative_chrX-5KB-sort.txt"]
    # npy_files = [f"KR_5kb_matrix_chr{chrom}_negative.npy" for chrom in range(1, 23)] + [
    #     "KR_5kb_matrix_chrX_negative.npy"]
    for txt_file, npy_file in zip(txt_files, npy_files):
        txt_path = os.path.join(txt_folder, txt_file)
        npy_path = os.path.join(npy_folder, npy_file)

        if os.path.exists(txt_path) and os.path.exists(npy_path):
            # Read the txt file to get chromosome data and coordinates
            chrom_data = []
            coord_data = []

            with open(txt_path, 'r') as f:
                for line in f:
                    columns = line.strip().split()
                    if len(columns) == 3:
                        chrom_data.append(columns[0])  # First column is chromosome as string
                        coord_data.append(
                            [float(columns[1]), float(columns[2])])  # Second and third columns are coordinates

            chrom_data = np.array(chrom_data).reshape(-1, 1)  # Convert to 2D array for concatenation
            coord_data = np.array(coord_data)

            # Read the .npy file
            npy_data = np.load(npy_path)

            # Combine chrom_data, coord_data with npy_data
            combined_data = np.hstack((chrom_data, coord_data, npy_data))

            # Define the output path
            output_path = os.path.join(output_folder, npy_file)

            # Save the combined data to the new output path
            np.save(output_path, combined_data)
            print(combined_data.shape)
            print(f"Updated and saved {npy_file} to {output_folder}")
        else:
            print(f"File not found: {txt_file} or {npy_file}")


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Merge txt and npy files and save them to a new directory.")

    # Command-line arguments
    parser.add_argument('-d', '--txt_folder', type=str, required=True,
                        help="Folder containing the txt files with coordinates.")
    parser.add_argument('-n', '--npy_folder', type=str, required=True, help="Folder containing the original npy files.")
    parser.add_argument('-o', '--output_folder', type=str, required=True,
                        help="Folder to save the merged npy files with coordinates.")

    # Parse arguments
    args = parser.parse_args()

    # Call the function with arguments from the command line
    add_txt_to_npy(args.txt_folder, args.npy_folder, args.output_folder)


if __name__ == '__main__':
    main()