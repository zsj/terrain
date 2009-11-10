/* ******************************************************************
 * Demo.cpp
 * 2009/11/07
 * zsj
 * *****************************************************************/

// INCLUDES /////////////////////////////////////////////////////////
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <GL/glut.h>

#pragma comment( linker, "/entry:\"mainCRTStartup\"" )

// DEFINES //////////////////////////////////////////////////////////
#define PI			3.14159265
#define EPSILON		0.000001

#define HMAP_SIZE	128				// (pixel) size of height map
#define BBOX_SIZE	16				// (pixel) size of bounding box

#define HMAP_INV	(1.0/(HMAP_SIZE-1))				// factor of size clamp to 1.0
#define VERTEX_NUM	(HMAP_SIZE*HMAP_SIZE)			// vertex number
#define FACE_NUM	((HMAP_SIZE-1)*(HMAP_SIZE-1)*2)	// surface number 

#define WALK_STEP	0.02			// walk step
#define MAX_NUM_COPLAN	6			// max number of coplanarity of one vertex

// MACROS ///////////////////////////
#define ATOR(a) ((a)/180.0*PI)		// converse from angle to radian
#define RTOA(r) ((r)/PI*180.0)		// converse form radian to angle

#define ADD(dest, v1, v2)		\
		dest.x = v1.x + v2.x;	\
		dest.y = v1.y + v2.y;	\
		dest.z = v1.z + v2.z;

#define SUB(dest, v1, v2)		\
		dest.x = v1.x - v2.x;	\
		dest.y = v1.y - v2.y;	\
		dest.z = v1.z - v2.z;

#define DIV(dest, v, factor)	\
		dest.x = v.x / factor;	\
		dest.y = v.y / factor;	\
		dest.z = v.z / factor;

#define MULT(dest, v, factor)	\
		dest.x = v.x * factor;	\
		dest.y = v.y * factor;	\
		dest.z = v.z * factor;

#define CROSS(dest, v1, v2)				\
		dest.x = v1.y*v2.z - v1.z*v2.y;	\
		dest.y = v1.z*v2.x - v1.x*v2.z;	\
		dest.z = v1.x*v2.y - v1.y*v2.x;

#define DOT(v1, v2)	(v1.x*v2.x+v1.y*v2.y+v1.z*v2.z)

#define NORMALIZE(n)	\
		nLenInv = 1.0 / sqrt(n.x*n.x + n.y*n.y + n.z*n.z);	\
		n.x *= nLenInv;	\
		n.y *= nLenInv;	\
		n.z *= nLenInv;

// STRUCTURES && CLASSES ////////////////////////////////////////////
struct VEC3I
{
	int x, y, z;
};

struct VEC3F
{
	float x, y, z;
};

struct SURFACE
{
	unsigned int vi[3];
};

struct BBOX
{
	VEC3F vmin, vmax;	// min, max inflexion coord of bounding box
	VEC3F center;		// center coord of bouning box
	float halflen[3];	// half length per axis
};

struct MESH
{
	unsigned int vnum;		// vertex number of mesh
	unsigned int fnum;		// surface number of mesh
	VEC3F		 *vlist;	// vertex list of mesh
	VEC3F		 *nvlist;	// normal list per vertex
	VEC3F		 *nflist;	// normal list per surface
	SURFACE		 *flist;	// surface list of mesh
	VEC3F		 center;	// center of mesh
};

struct HMAP						// height map structure
{
	float	*buffer;	// buffer of height map
	int		size;		// size of height map
	int		count;		// iterate count
	int		minDelta;	// min delta in height
	int		maxDelta;	// max delta in height
	float	filter;		// filter factor
	float	scaleX;		// scale factor of x aixs
	float	scaleY;		// scale factor of y aixs
	float	scaleZ;		// scale factor of z aixs
};

struct CAM						// camera structure
{
	VEC3F	vPosition;			// position 
	VEC3F	vLookat;			// look at direction
	VEC3F	vUp;				// up direction
	VEC3F	vForward;			// forward direction
	VEC3F	vSide;				// side direction
	float	vYaw;				// yaw
	float	vPitch;				// pitch
};

// GLOBALS //////////////////////////////////////////////////////////
float		  hmBuffer[HMAP_SIZE][HMAP_SIZE];		// height map buffer
unsigned char hmImage[HMAP_SIZE][HMAP_SIZE][3];		// height map image

VEC3F	vtList[VERTEX_NUM];		// vertex list
VEC3F	nvList[VERTEX_NUM];		// normal list per vertex
VEC3F	nfList[FACE_NUM];		// normal list per surface 
SURFACE	sfList[FACE_NUM];		// surface list

HMAP	htMap;					// height map
MESH	trnMesh;				// terrian mesh
CAM		camera;					// camera

int		mousePosX, mousePosY;	// mouse position

FILE	*log_file;				// log file
// FUNCTIONS ////////////////////////////////////////////////////////

// tools function /////////////////////////////////////////////
/* ***************************************************
 * print mesh
 * ***************************************************/
void printMesh(void)
{
	for( int i=0; i<VERTEX_NUM; i++ )
		fprintf(log_file, "%f %f %f\n", vtList[i].x, vtList[i].y, vtList[i].z);
}

/* ***************************************************
 * generate vector
 * ***************************************************/
VEC3F GenVector3f(float x, float y, float z)
{
	VEC3F vec;
	vec.x = x;
	vec.y = y;
	vec.z = z;
	return vec;
}

/////////////////////////////////////////////////////////////////////

/* ***************************************************
 * normalize height map
 * **************************************************/
void NormalizeImage(float *buffer, int size)
{
	int   bi;				// buffer index 
	float hValInv;			// interverse height value
	float minVal, maxVal;	// min and max value

	// init min && max value;
	minVal = buffer[0];
	maxVal = buffer[0];

	// find min && max value 
	for( bi=0; bi<size*size-1; bi++ )
	{
		if( buffer[bi]<minVal )
			minVal = buffer[bi];

		if( buffer[bi]>maxVal )
			maxVal = buffer[bi];
	}

	if( minVal>maxVal )
		return;

	hValInv = 1.0 / (maxVal - minVal);

	// normalize 
	for( bi=0; bi<size*size; bi++ )
		buffer[bi] = (buffer[bi]-minVal) * hValInv;
}

/* ***************************************************
 * filter band in height map
 * ***************************************************/
void FilterBand(float *band, int size, int stride, float factor)
{
	int bi, bj;			// band index 
	float prevVal;		// previous value 
	float factorBar;	// 1.0 - factor

	bj = stride;
	prevVal = band[0];
	factorBar = 1.0 - factor;

	for( bi=0; bi<size-1; bi++ )
	{
		band[bj] = prevVal*factor + band[bj]*factorBar;
		
		prevVal = band[bj];
		bj += stride;
	}
}

/* ***************************************************
 * filter height map
 * **************************************************/
void FilterImage(float *buffer, int size, float factor)
{
	int bi;		// buffer index
	
	// filter rows 
	for( bi=0; bi<size; bi++ )
	{
		// from left to right
		FilterBand(&buffer[bi*size], size, 1, factor);

		// from right to left
		FilterBand(&buffer[bi*size+size-1], size, -1, factor);
	}

	// filter columns
	for( bi=0; bi<size; bi++ )
	{
		// from top to bottom
		FilterBand(&buffer[bi], size, size, factor);

		// from bottom to top
		FilterBand(&buffer[size*(size-1)+bi], size, -size, factor);
	}
}

/* ***************************************************
 * generate height map using FAULT FORMATION ALGORITHM
 * in:	height map buffer pointer, width, height, 
 *		iterate count and begin and end increment 
 *		of height
 * ***************************************************/
void genHeightMap(HMAP *hMap)
{
	if( !hMap->buffer )
	{
		fprintf(log_file, "Invalidate height map! data is empty.");
		fclose(log_file);
		exit(0);
	}

	int bi;	// buffer index

	// initialize height map
	for( bi=0; bi<hMap->size*hMap->size; bi++ )
		hMap->buffer[bi] = 0.0;

	// iterate generate height map
	int   hDelta;		// increment of height
	float factor;		// factor of attenuation

	assert( hMap->count>0 );
	factor = (float)((hMap->minDelta-hMap->maxDelta) / hMap->count);

	int	ci;								// iterate count
	int	x, y;							// point in buffer
	unsigned int offset;				// position offset of buffer
	unsigned int x0, y0, x1, y1;		// begin and end point of line

	VEC3I	vl, vp;						// vector

	hDelta = hMap->maxDelta;			// init height increment
	
	for( ci=0; ci<hMap->count; ci++ )
	{
		// generate random line
		x0 = rand() % HMAP_SIZE;
		y0 = rand() % HMAP_SIZE;
		
		// obviate degenarate point
		do
		{
			x1 = rand() % HMAP_SIZE;
			y1 = rand() % HMAP_SIZE;
		} while( (x0==x1) && (y0==y1) );

		// vector of line
		vl.x = x1 - x0;
		vl.y = y1 - y0;
		vl.z = 0;

		// travel buffer and set point add dheight 
		// that lies left-side of line
		for( y=0; y<hMap->size; y++ )
		{
			vp.y   = y - y0;
			vp.z   = 0;
			offset = y*hMap->size;

			for( x=0; x<hMap->size; x++ )
			{			
				vp.x	= x - x0;
				if( (vl.x*vp.y - vl.y*vp.x)>0 )
					hMap->buffer[offset+x] += (float)hDelta;
			}
		}
		
		// filter 
		FilterImage(hMap->buffer, hMap->size, hMap->filter);

		// adjust height increment
		hDelta = hMap->maxDelta + ci*factor;
	}

	NormalizeImage(hMap->buffer, hMap->size);
}


/////////////////////////////////////////////////////////////////////

/* ******************************************************************
 * generate vertex info based on Height Map
 * *****************************************************************/
void genVertexList(HMAP *hMap)
{
	// vertex coord
	float fz;
	int x, z;
	unsigned int offset;	// position offset of buffer
	unsigned int vcount=0;	// vertex count
	
	// generate vertex coord
	for( z=0; z<hMap->size; z++ )
	{
		offset = z * hMap->size;
		fz = z * HMAP_INV * hMap->scaleZ;

		for( x=0; x<hMap->size; x++, vcount++ )
		{
			vtList[vcount].x = x * HMAP_INV * hMap->scaleX;
			vtList[vcount].y = hMap->buffer[offset+x] * hMap->scaleY;
			vtList[vcount].z = fz;
		}
	}
}

/* ******************************************************************
 * generate surface info based on Height Map
 * *****************************************************************/
void genSurfaceList(HMAP *hMap)
{
	// index of height map 
	int x, z;
	unsigned int offset;		// position offset of buffer
	unsigned int fcount = 0;	// surface count

	// generate surface info
	for( z=0; z<hMap->size-1; z++ )
	{
		offset = z * hMap->size;
		for( x=0; x<hMap->size-1; x++, fcount++ )
		{
			sfList[fcount].vi[0] = x + offset;
			sfList[fcount].vi[1] = x + offset + hMap->size + 1;
			sfList[fcount].vi[2] = x + offset + 1;

			fcount++;
			sfList[fcount].vi[0] = x + offset;
			sfList[fcount].vi[1] = x + offset + hMap->size;
			sfList[fcount].vi[2] = x + offset + hMap->size + 1;
		}
	}
}

/* ******************************************************************
 * generate normal(per surface) info
 * *****************************************************************/
void genSurfaceNormalList(void)
{
	unsigned int ni;	// index of normal list
	
	VEC3F ve1, ve2;		// vector of edge
	float nLenInv;		// inverse of normal length

	// generate normal per surface
	for( ni=0; ni<FACE_NUM; ni++ )
	{
		SUB(ve1, vtList[sfList[ni].vi[1]], vtList[sfList[ni].vi[0]]);
		SUB(ve2, vtList[sfList[ni].vi[2]], vtList[sfList[ni].vi[0]]);
		CROSS(nfList[ni], ve1, ve2);
		NORMALIZE(nfList[ni]);
	}
}

/* ******************************************************************
 * generate average normal of vertex that coplanarity surface
 * *****************************************************************/
VEC3F genAverageNormal(unsigned int *cpllist, unsigned int count)
{
	VEC3F vn;
	unsigned int ci;
	float nLenInv, countInv;

	countInv = 1.0 / count;
	vn.x = nfList[cpllist[0]].x;
	vn.y = nfList[cpllist[0]].y;
	vn.z = nfList[cpllist[0]].z;

	for( ci=1; ci<count; ci++ )
	{
		ADD(vn, vn, nfList[cpllist[ci]]); 
	}

	MULT(vn, vn, countInv);
	NORMALIZE(vn);

	return vn;
}

/* ******************************************************************
 * generate normal(per vertex) info 
 * *****************************************************************/
void genVertexNormalList(void)
{
	unsigned int ni = 0;							// index of normal list
	unsigned int rowi, coli;						// index of row and column of vertex grid
	unsigned int ci, ccount, cendi;					// index and count of coplanarity list
	unsigned int coplan_list[MAX_NUM_COPLAN];				// coplanarity surface list of one vertex
	unsigned int std_coplan_list[MAX_NUM_COPLAN*HMAP_SIZE];	// standard coplanarity list

	// first row ////////////////////////////////
	// first inflexion
	coplan_list[0]	= 0;
	coplan_list[1]	= 1;
	ccount			= 2;
	nvList[ni]		= genAverageNormal(coplan_list, ccount);
	
	// middle vertex normals of first row
	ni++;
	cendi	= ccount-1;
	ccount	= 3;
	for( coli=1; coli<HMAP_SIZE-1; coli++, ni++ )
	{
		coplan_list[0] = coplan_list[cendi];
		coplan_list[1] = coplan_list[cendi] + 1;
		coplan_list[2] = coplan_list[cendi] + 2;

		nvList[ni] = genAverageNormal(coplan_list, ccount);

		cendi = 2;
	}

	// final inflexion
	coplan_list[0]  = coplan_list[cendi];
	nvList[ni].x	= nfList[coplan_list[cendi]].x;
	nvList[ni].y	= nfList[coplan_list[cendi]].y;
	nvList[ni].z	= nfList[coplan_list[cendi]].z;

	// middle part ///////////////////////////////
	// second row
	ni++;
	coplan_list[0] = 0;
	coplan_list[1] = (HMAP_SIZE-1) * 2;
	coplan_list[2] = (HMAP_SIZE-1) * 2 + 1;
	ccount		   = 3;
	nvList[ni]	   = genAverageNormal(coplan_list, ccount);

	// copy first column coplanarity list to standard coplan list
	memcpy((void*)std_coplan_list, (void*)coplan_list, MAX_NUM_COPLAN*sizeof(unsigned int));

	// middle vertex normals of second row
	ni++;
	cendi	= ccount-1;
	ccount	= MAX_NUM_COPLAN;
	for( coli=1, ci=MAX_NUM_COPLAN; coli<HMAP_SIZE-1; coli++, ni++, ci+=MAX_NUM_COPLAN )
	{
		coplan_list[0] = coplan_list[0];
		coplan_list[1] = coplan_list[0] + 1;
		coplan_list[2] = coplan_list[0] + 2;
		coplan_list[3] = coplan_list[cendi];
		coplan_list[4] = coplan_list[cendi] + 1;
		coplan_list[5] = coplan_list[cendi] + 2;

		nvList[ni] = genAverageNormal(coplan_list, ccount);

		// copy current column coplan list to standard coplan list
		memcpy((void*)&std_coplan_list[ci], (void*)coplan_list, MAX_NUM_COPLAN*sizeof(unsigned int));

		cendi = 5;
	}

	// final vertex's normal of second row
	coplan_list[0] = (HMAP_SIZE-1) * 2 - 2;
	coplan_list[1] = (HMAP_SIZE-1) * 2 - 1;
	coplan_list[2] = (HMAP_SIZE-1) * 4 - 1;
	ccount		   = 3;
	nvList[ni]	   = genAverageNormal(coplan_list, ccount);

	// copy final column coplan list to standard coplan list
	memcpy((void*)&std_coplan_list[ci-MAX_NUM_COPLAN], (void*)coplan_list, MAX_NUM_COPLAN*sizeof(unsigned int));

	// other rows in middle part
	ni++;
	for( rowi=2; rowi<HMAP_SIZE-1; rowi++ )
	{
		// adjust standard coplan list
		for( ci=0; ci<MAX_NUM_COPLAN*HMAP_SIZE; ci++ )
			std_coplan_list[ci] += (HMAP_SIZE-1)*2;

		// first vertex's normal of current row
		ccount = 3;
		nvList[ni] = genAverageNormal(std_coplan_list, ccount);

		// middle vertices normal of current row
		ni++;
		ccount = MAX_NUM_COPLAN;
		for( coli=1, ci=MAX_NUM_COPLAN; coli<HMAP_SIZE-1; coli++, ni++, ci+=MAX_NUM_COPLAN )
			nvList[ni] = genAverageNormal(&std_coplan_list[ci], ccount);

		// final vertex's normal of current row
		ccount = 3;
		nvList[ni] = genAverageNormal(&std_coplan_list[ci-MAX_NUM_COPLAN], ccount);

		// adjust standard coplan list
		for( ci=3; ci<MAX_NUM_COPLAN; ci++ )
			std_coplan_list[ci] = 0;
		for( ci=MAX_NUM_COPLAN*HMAP_SIZE-1; ci>=MAX_NUM_COPLAN*HMAP_SIZE-3; ci-- )
			std_coplan_list[ci] = 0;
	}

	// final row ////////////////////////////////
	// first vertex's normal
	ni++;
	nvList[ni].x = nfList[(HMAP_SIZE-1) * (HMAP_SIZE-1) * 2].x;
	nvList[ni].y = nfList[(HMAP_SIZE-1) * (HMAP_SIZE-1) * 2].y;
	nvList[ni].z = nfList[(HMAP_SIZE-1) * (HMAP_SIZE-1) * 2].z;

	// other vertices normal of middle part in final row
	ni++;
	cendi  = 0;
	ccount = 3;
	for( coli=1; coli<HMAP_SIZE-1; coli++, ni++ )
	{
		coplan_list[0] = coplan_list[cendi];
		coplan_list[1] = coplan_list[cendi] + 1;
		coplan_list[2] = coplan_list[cendi] + 2;

		nvList[ni] = genAverageNormal(coplan_list, ccount);

		cendi = 2;
	}

	// final vertex's normal
	coplan_list[0] = coplan_list[cendi];
	coplan_list[1] = coplan_list[cendi] + 1;
	ccount		   = 2;
	nvList[ni]	   = genAverageNormal(coplan_list, ccount);
}

/////////////////////////////////////////////////////////////////////

// setup camera
void SetupCamera(void)
{
	float cosYaw = (float)cos(ATOR(camera.vYaw));
	float sinYaw = (float)sin(ATOR(camera.vYaw));
	float cosPitch = (float)cos(ATOR(camera.vPitch));
	float sinPitch = (float)sin(ATOR(camera.vPitch));

	// calculate forward direction
	camera.vForward.x = sinYaw * cosPitch;
	camera.vForward.y = sinPitch;
	camera.vForward.z = cosPitch * -cosYaw;

	// calculate lookat direction
	ADD(camera.vLookat, camera.vPosition, camera.vForward);

	// calculate side direction
	CROSS(camera.vSide, camera.vForward, camera.vUp);

	// setup camera
	gluLookAt(camera.vPosition.x, camera.vPosition.y, camera.vPosition.z,
			  camera.vLookat.x, camera.vLookat.y, camera.vLookat.z,
			  camera.vUp.x, camera.vUp.y, camera.vUp.z);
}

// init
void init(void)
{	
	if( !(log_file = fopen("log.txt", "wb")) )
		exit(0);
	
	//srand((unsigned int)time(NULL));

	// initialize height map
	htMap.buffer	= &hmBuffer[0][0];
	htMap.size		= HMAP_SIZE;
	htMap.count		= 64;
	htMap.maxDelta	= 128;
	htMap.minDelta	= 0;
	htMap.filter	= 0.2;
	htMap.scaleX	= 4.0;
	htMap.scaleY	= 1.0;
	htMap.scaleZ	= 4.0;

	// generate height map
	genHeightMap(&htMap);
	//filterImage(0.8);

	// initialize terrian mesh
	trnMesh.vnum	= VERTEX_NUM;
	trnMesh.fnum	= FACE_NUM;
	trnMesh.vlist	= vtList;
	trnMesh.flist	= sfList;
	trnMesh.nvlist  = nvList;
	trnMesh.nflist  = nfList;
	trnMesh.center.x = 0.5;
	trnMesh.center.y = 0;
	trnMesh.center.z = 0.5;

	// generate terrian mesh
	genVertexList(&htMap);
	genSurfaceList(&htMap);
	genSurfaceNormalList();
	genVertexNormalList();

	// initialize camera	
	camera.vPosition = GenVector3f(1.0, 1.0, 3.5);
	camera.vUp		 = GenVector3f(0.0, 1.0, 0.0);
	camera.vSide	 = GenVector3f(1.0, 0.0, 0.0);
	camera.vYaw		 = 0.0f;
	camera.vPitch	 = 0.0f;

	// set lighting information
	GLfloat mat_specular[]		= {0.0, 1.0, 0.0, 1.0};
	GLfloat mat_diffuse[]		= {0.0, 0.8, 0.0, 0.5};
	GLfloat mat_shininess[]		= {50.0};
	GLfloat light_position[]	= {1.0, 1.0, 1.0, 0.0};
	GLfloat white_light[]		= {1.0, 1.0, 1.0, 1.0};
	GLfloat diffuse_light[]		= {0.8, 0.8, 0.8, 0.8};
	GLfloat lmodel_ambient[]	= {0.8, 1.0, 0.8, 0.8};

	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse_light);
	glLightfv(GL_LIGHT0, GL_SPECULAR, white_light);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);

	// enable lighting
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);

	// close log file
	fclose(log_file);
}

// display
void display(void)
{
	glClearColor(0.18, 0.18, 0.776, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_CULL_FACE);	
	glCullFace(GL_BACK);
	glPolygonMode(GL_FRONT, GL_FILL);	

	glPushMatrix();
	glLoadIdentity();

	// setup camera
	SetupCamera();

	glShadeModel(GL_SMOOTH);
	glColor3f(0.0, 0.8, 0.0);
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	
	glVertexPointer(3, GL_FLOAT, 0, vtList);
	glNormalPointer(GL_FLOAT, 0, nvList);

	glDrawElements(GL_TRIANGLES, FACE_NUM*3, GL_UNSIGNED_INT, sfList);

	glPopMatrix();

	glutSwapBuffers();
}

// reshape
void reshape(int w, int h)
{
	if( h<1 )
		h=1;

	glViewport(0, 0, (GLsizei) w, (GLsizei) h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(90.0, (GLfloat)(w)/(GLfloat)(h-21), 0.1, 100.0);
	glMatrixMode(GL_MODELVIEW);
}

// keyboard
void keyboard(unsigned char key, int x, int y)
{
	switch(key)
	{
	case 'w':
		VEC3F vec;
		MULT(vec, camera.vForward, WALK_STEP);
		ADD(camera.vPosition, camera.vPosition, vec);
		glutPostRedisplay();
		break;

	case 's':
		MULT(vec, camera.vForward, WALK_STEP);
		SUB(camera.vPosition, camera.vPosition, vec);
		glutPostRedisplay();
		break;

	case 'd':
		MULT(vec, camera.vSide, WALK_STEP);
		ADD(camera.vPosition, camera.vPosition, vec);
		glutPostRedisplay();
		break;

	case 'a':
		MULT(vec, camera.vSide, WALK_STEP);
		SUB(camera.vPosition, camera.vPosition, vec);
		glutPostRedisplay();
		break;

	case 'q':
		// free memory

		exit(0);
		break;
	}
}

// mouse func
void mouse(int button, int state, int x, int y)
{
	mousePosX = x;
	mousePosY = y;
}

// motion func
void motion(int x, int y)
{
	camera.vYaw	  += (float)(x-mousePosX) * 0.1f;
	camera.vPitch -= (float)(y-mousePosY) * 0.1f;

	if( camera.vYaw>=360.0 || camera.vYaw<=-360.0 )
		camera.vYaw=0.0;

	if( camera.vPitch>60.0 )
		camera.vPitch=60.0;

	if( camera.vPitch<-60.0 )
		camera.vPitch=-60.0;

	glutPostRedisplay();
}

// ENTRY ////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(800, 600);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("Terrian Demo");
	init();
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	glutMainLoop();
	return 0;
}