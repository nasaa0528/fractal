// Fractal.cpp : Defines the entry point for the console application.
//

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <complex>
#include <mpi.h>
using namespace std;

void WriteTGA_RGB(const char* filename, unsigned char* data, unsigned int width, unsigned int height)
{
	FILE *f = fopen(filename, "wb");
	if (!f) {
		fprintf(stderr, "Unable to create output TGA image `%s'\n", filename);
		exit(EXIT_FAILURE);
	}

	fputc(0x00, f); /* ID Length, 0 => No ID        */
	fputc(0x00, f); /* Color Map Type, 0 => No color map included   */
	fputc(0x02, f); /* Image Type, 2 => Uncompressed, True-color Image */
	fputc(0x00, f); /* Next five bytes are about the color map entries */
	fputc(0x00, f); /* 2 bytes Index, 2 bytes length, 1 byte size */
	fputc(0x00, f);
	fputc(0x00, f);
	fputc(0x00, f);
	fputc(0x00, f); /* X-origin of Image    */
	fputc(0x00, f);
	fputc(0x00, f); /* Y-origin of Image    */
	fputc(0x00, f);
	fputc(width & 0xff, f); /* Image Width      */
	fputc((width >> 8) & 0xff, f);
	fputc(height & 0xff, f); /* Image Height     */
	fputc((height >> 8) & 0xff, f);
	fputc(0x18, f); /* Pixel Depth, 0x18 => 24 Bits */
	fputc(0x20, f); /* Image Descriptor     */

	for (int y = height - 1; y >= 0; y--) {
		for (size_t x = 0; x < width; x++) {
			const size_t i = (y * width + x) * 3;
			fputc(data[i + 2], f); /* write blue */
			fputc(data[i + 1], f); /* write green */
			fputc(data[i], f); /* write red */
		}
	}
}

int main(int argc, char **argv)
{
	// MPI declarations 
	int id, nproc; 
	MPI_Status status; 

	const unsigned int domainWidth = 1024;
	const unsigned int domainHeight = 1024;

	complex<double> K(0.353, 0.288);
	complex<double> center(-1.68, -1.23);
	double scale = 2.35;
	const unsigned int maxIterations = 100;


	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	unsigned char *data = new unsigned char[domainWidth * domainHeight * 3];
	memset(data, 0, domainWidth * domainHeight * 3 * sizeof(char));

	if(id == 0){
		// cout << "ID is:\t" << id << "\tsize:\t" << nproc << "\tdomainWidth:\t" << domainWidth << "\tdomainHeight:\t" << domainHeight << "\tscale:\t" << scale << "\tmaxIterations:\t" << maxIterations << "\tK:\t" << K << "\tcenter:\t" << center << "\n\n" << endl; 
		
	/************* MASTER PROCESS DOES ITS JOB STARTS HERE ******************/
		// for (unsigned int y = 0; y < domainHeight; ++y)
		for (unsigned int y = id; y < domainHeight; y+=nproc)
		{
			for (unsigned int x = 0; x < domainWidth; ++x)
			{
				complex<double> c(x / (double)domainWidth * scale + center.real(),
					y / (double)domainHeight * scale + center.imag());

				complex<double> z(c);
				for (unsigned int iteration = 0; iteration < maxIterations; ++iteration)
				{
					z = z * z + c;
					if (abs(z) > 1.0f)
					{
						data[(x + y * domainWidth) * 3 + 0] = 255;
						data[(x + y * domainWidth) * 3 + 1] = 255;
						data[(x + y * domainWidth) * 3 + 2] = 255;
					}
				}
			}
		}
	/************* MASTER PROCESS DOES ITS JOB ENDS HERE ******************/

	/************* MASTER RECEIVE AND FILL CONCATENATE JOB STARTS HERE ******************/
		for (unsigned int rank_id = 1; rank_id < nproc; rank_id++){
			unsigned char *rec_buf = new unsigned char[domainWidth * domainHeight * 3];
			memset(rec_buf, 0, domainWidth * domainHeight * 3 * sizeof(char));
			MPI_Recv(rec_buf, (domainWidth * domainHeight * 3), MPI_UNSIGNED_CHAR, rank_id, 1, MPI_COMM_WORLD, &status);


			for (unsigned int y = rank_id; y < domainHeight; y+=nproc)
			{
				for (unsigned int x = 0; x < domainWidth; ++x)
				{
					data[(x + y * domainWidth) * 3 + 0] = rec_buf[(x + y * domainWidth) * 3 + 0];
					data[(x + y * domainWidth) * 3 + 1] = rec_buf[(x + y * domainWidth) * 3 + 1];
					data[(x + y * domainWidth) * 3 + 2] = rec_buf[(x + y * domainWidth) * 3 + 2];
				}
			}

			delete[] rec_buf;
		}

	/************* MASTER RECEIVE AND FILL CONCATENATE JOB ENDS HERE ******************/

		WriteTGA_RGB("mandelbrot_par.tga", data, domainWidth, domainHeight);
		delete[] data;
	}
	else { 
		// cout << "ID is:\t" << id << "\tsize:\t" << nproc << "\tdomainWidth:\t" << domainWidth << "\tdomainHeight:\t" << domainHeight << "\tscale:\t" << scale << "\tmaxIterations:\t" << maxIterations << "\tK:\t" << K << "\tcenter:\t" << center << "\n\n" << endl; 

		/************* SLAVE PROCESS DOES ITS JOB STARTS HERE ******************/
		// for (unsigned int y = 0; y < domainHeight; ++y)
		for (unsigned int y = id; y < domainHeight; y+=nproc)
		{
			for (unsigned int x = 0; x < domainWidth; ++x)
			{
				complex<double> c(x / (double)domainWidth * scale + center.real(),
					y / (double)domainHeight * scale + center.imag());

				complex<double> z(c);
				for (unsigned int iteration = 0; iteration < maxIterations; ++iteration)
				{
					z = z * z + c;
					if (abs(z) > 1.0f)
					{
						data[(x + y * domainWidth) * 3 + 0] = 255;
						data[(x + y * domainWidth) * 3 + 1] = 255;
						data[(x + y * domainWidth) * 3 + 2] = 255;
					}
				}
			}
		}
		MPI_Send(data, (domainWidth * domainHeight * 3), MPI_UNSIGNED_CHAR, 0, 1, MPI_COMM_WORLD);

		/************* SLAVE PROCESS DOES ITS JOB ENDS HERE ******************/

		// cout << "if the data and tmp are equal" << strcmp(tmp,data) << endl;
	}
	MPI_Finalize();
	return 0;
}

