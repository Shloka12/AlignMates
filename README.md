# AlignMates
AlignMates is an interactive application designed to visualize the Needleman-Wunsch algorithm, a fundamental algorithm used in bioinformatics for aligning genomic/proteomic sequences. It allows users to input two sequences and set parameters for match, mismatch, and gap penalties to compute the optimal global alignment.

<img width="749" alt="Welcome" src="https://github.com/Shloka12/AlignMates/assets/67782856/845cc28d-1f2b-410a-a56f-6094752a3ed5">

<img width="745" alt="help" src="https://github.com/Shloka12/AlignMates/assets/67782856/673e738a-8160-4238-869e-49e539fcf6b5">

<img width="746" alt="Alignment" src="https://github.com/Shloka12/AlignMates/assets/67782856/87b94a97-0fcf-474e-92ca-15e83a2094fe">


### Prerequisites
- Python 3.x
- cmu_graphics
- PIL (Python Imaging Library)

## Project Structure
- `src/`: Contains the Python script (`AlignMates.py`) for running the application.
- `images/`: Contains image files (`dna.png`, `help.png`) used by the application.
- `lib/`: Contains any libraries or additional code required for the application.

## Usage
- **Tab**: Switches between input fields.
- **Enter**: Executes the sequence alignment algorithm.
- **Escape/Space/Right Arrow**: Switches between screens/views.

## Additional Notes
Ensure all input sequences are valid DNA/RNA/protein sequences. The application is case-insensitive for sequence inputs. Match, mismatch, and gap scores should be integers.


