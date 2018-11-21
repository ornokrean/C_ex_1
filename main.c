/**
 * @file AnalyzeProtein.c
 * @author  Or Nokrean <ornokrean@mail.huji.ac.il>
 * @version 1.0
 * @date 31 Oct 2018
 *
 * @brief count atoms from given files and analyze them.
 *
 * @section LICENSE
 * This program is not a free software; contact Or Nokrean for license.
 *
 * @section DESCRIPTION
 * The program analyzes the given file and outputs the atoms data.
 * Input  : atom file(s) (PDB).
 * Process: counting the atoms, saving them and calculating the center gravity and other
 *          element from the data.
 * Output : printed data that was analyzed from the given file(s).
 */

// ------------------------------ includes ------------------------------
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
// -------------------------- const definitions -------------------------

/**
 * @def X_INDEX 0
 * @brief The index of the x coordinate in each atom
 */
#define X_INDEX 0

/**
 * @def Y_INDEX 1
 * @brief The index of the y coordinate in each atom
 */
#define Y_INDEX 1

/**
 * @def Z_INDEX 2
 * @brief The index of the z coordinate in each atom
 */
#define Z_INDEX 2

/**
 * @def X_STRING_INDEX 30
 * @brief The index of the x coordinate substring in each atom line
 */
#define X_STRING_INDEX 30

/**
 * @def Y_STRING_INDEX 38
 * @brief The index of the y coordinate substring in each atom line
 */
#define Y_STRING_INDEX 38

/**
 * @def Z_STRING_INDEX 46
 * @brief The index of the z coordinate substring in each atom line
 */
#define Z_STRING_INDEX 46

/**
 * @def NUM_OF_CORD 3
 * @brief The number of coordinates in each atom
 */
#define NUM_OF_CORD 3

/**
 * @def COORDINATE_CHAR_LEN 8
 * @brief The number of coordinates in each atom
 */
#define COORDINATE_CHAR_LEN 8

/**
 * @def LINE_MIN_LEN 61
 * @brief The minimum atom line length
 */
#define LINE_MIN_LEN 61

/**
 * @def LINE_MAX_LEN 80
 * @brief The maximum atom line length
 */
#define LINE_MAX_LEN 80

/**
 * @def MAX_NUMBER_OF_ATOMS 20000
 * @brief The maximum number of atoms in a file
 */
#define MAX_NUMBER_OF_ATOMS 20000

/**
 * @def MIN_ARGS_NUM 1
 * @brief The minimum number of args in a program run
 */
#define MIN_ARGS_NUM 1


/**
 * @var char ATOM[]
 * @brief Start of a valid atom line
 */
const char ATOM[] = "ATOM  ";

/**
 * @var char COORDINATE_CONVERSION_ERR[]
 * @brief Error in coordinate conversion message
 */
const char COORDINATE_CONVERSION_ERR[] = "Error in coordinate conversion %s!\n";

/**
 * @var char SHORT_ATOM_LINE_ERR[]
 * @brief Atom line too short error message
 */
const char SHORT_ATOM_LINE_ERR[] = "ATOM line is too short %ld characters\n";

/**
 * @var char OPENING_FILE_ERR[]
 * @brief Error in opening file message
 */
const char OPENING_FILE_ERR[] = "Error opening file: %s\n";

/**
 * @var char ZERO_ATOMS_ERR[]
 * @brief No valid atom(s) message
 */
const char ZERO_ATOMS_ERR[] = "Error - 0 atoms were found in the file %s\n";

/**
 * @var char WRONG_USAGE_ERR[]
 * @brief No valid atom(s) message
 */
const char WRONG_USAGE_ERR[] = "Usage: AnalyzeProtein <pdb1> <pdb2> ...\n";

/**
 * @var char PDB_FILE_NAME_AND_ATOMS_OUT[]
 * @brief the first line of the results, the file name and the number of valid atoms
 *          found in the file
 */
const char PDB_FILE_NAME_AND_ATOMS_OUT[] = "PDB file %s, %d atoms were read\n";

/**
 * @var char CG_OUT[]
 * @brief the second line of the results, the center gravity of the found atoms.
 */
const char CG_OUT[] = "Cg = %.3f %.3f %.3f\n";

/**
 * @var char RG_OUT[]
 * @brief the third line of the results, the Rg of the found atoms.
 */
const char RG_OUT[] = "Rg = %.3f\n";

/**
 * @var char DMAX_OUT[]
 * @brief the 4th line of the results, the Dmax of the found atoms.
 */
const char DMAX_OUT[] = "Dmax = %.3f\n";


// ------------------------------ functions -----------------------------

/**
 * @brief function to calculate the distance between two given points
 * @param point_one the first point
 * @param point_two the second point
 * @return float - the distance between the two points
 */
float getDistance(float point_one[NUM_OF_CORD], float point_two[NUM_OF_CORD])
{
    return sqrtf(powf(point_one[X_INDEX] - point_two[X_INDEX], 2) +
                 powf(point_one[Y_INDEX] - point_two[Y_INDEX], 2) +
                 powf(point_one[Z_INDEX] - point_two[Z_INDEX], 2));
}

/**
 * @brief function to calculate the Center Gravity of the atoms given
 * @param allAtoms all of the atoms read from the file
 * @param numberOfAtoms the number of atoms in allAtoms
 * @param result the result of the center gravity, written into the given pointer to array
 */
void getCenterGravity(float allAtoms[][NUM_OF_CORD], int numberOfAtoms, float *result)
{
    result[X_INDEX] = 0;
    result[Y_INDEX] = 0;
    result[Z_INDEX] = 0;
    for (int i = 0; i < numberOfAtoms; i++)
    {
        result[X_INDEX] += allAtoms[i][X_INDEX];
        result[Y_INDEX] += allAtoms[i][Y_INDEX];
        result[Z_INDEX] += allAtoms[i][Z_INDEX];
    }
    result[X_INDEX] /= numberOfAtoms;
    result[Y_INDEX] /= numberOfAtoms;
    result[Z_INDEX] /= numberOfAtoms;
}

/**
 * @brief this function get a certain coordinate from the given line by the given boundaries
 *            from index to COORDINATE_CHAR_LEN
 * @param line the line to take coordinates from
 * @param index the start index of the coordinate to take in the line
 * @return the float value of the wanted coordinate, or exit if the float conversion had a problem
 */
float getCoordinate(char *line, int index)
{
    char coordinateSubstr[COORDINATE_CHAR_LEN], *end = NULL, *actualCoordinate = NULL;
    float numericCoordinate = 0;

    // Get the substring of the coordinate
    actualCoordinate = strncpy(coordinateSubstr, &line[index], COORDINATE_CHAR_LEN);

    // make it a float and check for valid float
    numericCoordinate = strtof(actualCoordinate, &end);
    errno = 0;
    if (numericCoordinate == 0 && (errno != 0 || end == actualCoordinate))
    {
        fprintf(stderr, COORDINATE_CONVERSION_ERR, actualCoordinate);
        exit(EXIT_FAILURE);
    }
    return numericCoordinate;
}

/**
 * @brief this function adds an atom to allAtoms list if the line starts with "ATOM" and is
 *          between 60 to 80 chars
 * @param allAtoms all of the atom that was already read from the file
 * @param numberOfAtoms the number of current atom that were read from the file
 * @param line the line to check if it contains atom and if so, to read and add the coordinates
 *          using getCoordinate
 * @return the number of atoms after checking if it is a real atom and adding it to the list
 */
int addAtom(float allAtoms[][NUM_OF_CORD], int numberOfAtoms, char line[LINE_MAX_LEN])
{
    // line starts with "ATOM"
    if (strstr(line, ATOM) == line)
    {
        // line is no shorter then 60 chars
        size_t lineLength = strlen(line);
        if (lineLength < LINE_MIN_LEN)
        {
            fprintf(stderr, SHORT_ATOM_LINE_ERR, lineLength);
            exit(EXIT_FAILURE);
        }
        //add atom cords
        allAtoms[numberOfAtoms][X_INDEX] = getCoordinate(line, X_STRING_INDEX);
        allAtoms[numberOfAtoms][Y_INDEX] = getCoordinate(line, Y_STRING_INDEX);
        allAtoms[numberOfAtoms][Z_INDEX] = getCoordinate(line, Z_STRING_INDEX);
        numberOfAtoms++;
    }
    return numberOfAtoms;

}

/**
 * @brief this function gets the data from the file, by reading it line by line and adding the
 *          atom if such were found
 * @param allAtoms all the atoms read from the file
 * @param filePath the address of the file to read from
 * @return the number of atoms that were correctly read from the file
 */
int getFileData(float allAtoms[][NUM_OF_CORD], char *filePath)
{
    FILE *file;
    file = fopen(filePath, "r");
    if (file == NULL)
    {
        fprintf(stderr, OPENING_FILE_ERR, filePath);
        exit(EXIT_FAILURE);
    }
    int numberOfAtoms = 0;
    char line[LINE_MAX_LEN];
    while (fgets(line, LINE_MAX_LEN, file) != NULL)
    {
        // check for more than 20,000 atoms
        if (numberOfAtoms >= MAX_NUMBER_OF_ATOMS)
        {
            break;
        }
        numberOfAtoms = addAtom(allAtoms, numberOfAtoms, line);
    }
    fclose(file);
    return numberOfAtoms;
}


/**
 * @brief this function prints the results of the file that were given in filePath
 * @param allAtoms all the atoms that were read from the file
 * @param numberOfAtoms the number of correctly read atoms from the file
 * @param filePath the address of the wanted file
 */
void printFileData(float allAtoms[][NUM_OF_CORD], int numberOfAtoms, char *filePath)
{
    // get the center gravity
    float centerGravity[NUM_OF_CORD];
    getCenterGravity(allAtoms, numberOfAtoms, centerGravity);

    //get the rg and the dmax
    float rg = 0, dmax = 0, tempDistance = 0;
    for (int i = 0; i < numberOfAtoms; i++)
    {
        for (int j = i + 1; j < numberOfAtoms; j++)
        {
            tempDistance = getDistance(allAtoms[i], allAtoms[j]);
            if (tempDistance > dmax)
            {
                dmax = tempDistance;
            }

        }
        rg += powf(getDistance(centerGravity, allAtoms[i]), 2);
    }
    rg = sqrtf(rg / numberOfAtoms);

    // print all the data
    printf(PDB_FILE_NAME_AND_ATOMS_OUT, filePath, numberOfAtoms);
    printf(CG_OUT, centerGravity[X_INDEX], centerGravity[Y_INDEX], centerGravity[Z_INDEX]);
    printf(RG_OUT, rg);
    printf(DMAX_OUT, dmax);
}


/**
 * @brief this function runs the getFileData and checks if any atoms were read, if so, it will
 *          be printed by printFileData function, else it will print a warning
 * @param filePath the file to read from
 */
void fileAnalyze(char *filePath)
{
    float allAtoms[MAX_NUMBER_OF_ATOMS][NUM_OF_CORD];
    int numberOfAtoms = getFileData(allAtoms, filePath);
    if (numberOfAtoms == 0)
    {
        fprintf(stderr, ZERO_ATOMS_ERR, filePath);
        exit(EXIT_FAILURE);
    }

    printFileData(allAtoms, numberOfAtoms, filePath);

}

/**
 * @brief this is the main function of the program
 * @param argc the number f args given
 * @param argv the array contains all the args given
 * @return 1 for error, 0 for success.
 */
int main(int argc, char *argv[])
{
    if (argc <= MIN_ARGS_NUM)
    {
        fprintf(stderr, WRONG_USAGE_ERR);
        return 1;
    }
    for (int i = 1; i < argc; i++)
    {
        fileAnalyze(argv[i]);
    }
    return 0;
}

