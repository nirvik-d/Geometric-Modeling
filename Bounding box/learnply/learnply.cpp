/*

Functions for learnply

Eugene Zhang, 2005
*/
#include<windows.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "glut.h"
#include <string>
#include <fstream>
# include <sstream>
#include <iostream>
#include "ply.h"
#include "icVector.H"
#include "icMatrix.H"
#include "learnply_io.h"
#include "trackball.h"
#include "tmatrix.h"
#include "jacobi.cpp"
#include "nrutil.cpp"
#include "learnply.h"

//#define NRANSI

using namespace std;

static PlyFile *in_ply;

unsigned char orientation;  // 0=ccw, 1=cw

FILE *this_file;
const int win_width = 1024;
const int win_height = 1024;

double radius_factor = 0.9;

int display_mode = 0;
double error_threshold = 1.0e-13;
char reg_model_name[128];
FILE *f;
int ACSIZE = 1; // for antialiasing
int view_mode = 0;  // 0 = othogonal, 1=perspective
float s_old, t_old;
float rotmat[4][4];
static Quaternion rvec;

int mouse_mode = -2;  // -2=no action, -1 = down, 0 = zoom, 1 = rotate x, 2 = rotate y, 3 = tranlate x, 4 = translate y, 5 = cull near 6 = cull far
int mouse_button = -1; // -1=no button, 0=left, 1=middle, 2=right
int last_x, last_y;

struct jitter_struct {
	double x;
	double y;
} jitter_para;

jitter_struct ji1[1] = { { 0.0, 0.0 } };
jitter_struct ji16[16] = { { 0.125, 0.125 },{ 0.375, 0.125 },{ 0.625, 0.125 },{ 0.875, 0.125 },
{ 0.125, 0.375 },{ 0.375, 0.375 },{ 0.625, 0.375 },{ 0.875, 0.375 },
{ 0.125, 0.625 },{ 0.375, 0.625 },{ 0.625, 0.625 },{ 0.875, 0.625 },
{ 0.125, 0.875 },{ 0.375, 0.875 },{ 0.625, 0.875 },{ 0.875, 0.875 }, };

Polyhedron *poly;

void init(void);
void keyboard(unsigned char key, int x, int y);
void motion(int x, int y);
void display(void);
void mouse(int button, int state, int x, int y);
void display_shape(GLenum mode, Polyhedron *poly);


void axisBoundingBox(Polyhedron *this_poly);
void momentBoundingBox(Polyhedron *this_poly);
void normBoundingBox(Polyhedron *this_poly);
void drawLine(icVector3 &x, icVector3 &y, float width, float r, float g, float b);
void drawBoundingBox(Polyhedron *this_poly);


/******************************************************************************
Main program.
******************************************************************************/
int operationNum = -1; //tell program to do which operation
float L = 1.0;//parameter L of checkerboard
bool Lflag = true;

int main(int argc, char *argv[])
{
	char *progname;
	int num = 1;
	FILE *this_file;

	progname = argv[0];

	this_file = fopen("../tempmodels/bunny.ply", "r");
	poly = new Polyhedron(this_file);
	fclose(this_file);
	mat_ident(rotmat);

	poly->initialize(); // initialize everything

	poly->calc_bounding_sphere();
	poly->calc_face_normals_and_area();
	poly->average_normals();

	cout << "Press 1 to generate axis aligned bounding box" << endl;
	cout << "Press 2 to generate moment based bounding box" << endl;
	cout << "Press 3 to generate norm based bounding box" << endl;

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowPosition(20, 20);
	glutInitWindowSize(win_width, win_height);
	glutCreateWindow("Geometric Modeling");
	init();
	glutKeyboardFunc(keyboard);
	glutDisplayFunc(display);
	glutMotionFunc(motion);
	glutMouseFunc(mouse);
	glutMainLoop();
	poly->finalize();  // finalize everything

					   //}
					   /***********************************************************************/


	return 0;    /* ANSI C requires main to return int. */
}

void color_mapping(double percentage, double col[3])
{
	if (percentage == 0.0) {
		col[0] = 1.0;
		col[1] = 1.0;
		col[2] = 1.0;
	}
	else if (percentage <= 1.0 / 3) {
		col[0] = 1.0;
		col[1] = 1.0 - percentage*3.0;
		col[2] = 1.0 - percentage*3.0;
	}
	else if (percentage <= 2.0 / 3) {
		col[0] = 1.0;
		col[1] = percentage*3.0 - 1.0;
		col[2] = 0.0;
	}
	else if (percentage <= 3.0 / 3) {
		col[0] = 3.0 - percentage*3.0;
		col[1] = 1.0;
		col[2] = 0.0;
	}
	else {
		col[0] = 1.0;
		col[1] = 1.0;
		col[2] = 0.0;
	}
}

/******************************************************************************
Read in a polyhedron from a file.
******************************************************************************/

Polyhedron::Polyhedron(FILE *file)
{
	int i, j;
	int elem_count;
	char *elem_name;

	/*** Read in the original PLY object ***/
	in_ply = read_ply(file);

	for (i = 0; i < in_ply->num_elem_types; i++) {

		/* prepare to read the i'th list of elements */
		elem_name = setup_element_read_ply(in_ply, i, &elem_count);

		if (equal_strings("vertex", elem_name)) {

			/* create a vertex list to hold all the vertices */
			nverts = max_verts = elem_count;
			vlist = new Vertex *[nverts];

			/* set up for getting vertex elements */

			setup_property_ply(in_ply, &vert_props[0]);
			setup_property_ply(in_ply, &vert_props[1]);
			setup_property_ply(in_ply, &vert_props[2]);
			vert_other = get_other_properties_ply(in_ply,
				offsetof(Vertex_io, other_props));

			/* grab all the vertex elements */
			for (j = 0; j < nverts; j++) {
				Vertex_io vert;
				get_element_ply(in_ply, (void *)&vert);

				/* copy info from the "vert" structure */
				vlist[j] = new Vertex(vert.x, vert.y, vert.z);
				vlist[j]->other_props = vert.other_props;
			}
		}
		else if (equal_strings("face", elem_name)) {

			/* create a list to hold all the face elements */
			ntris = max_tris = elem_count;
			tlist = new Triangle *[ntris];

			/* set up for getting face elements */
			setup_property_ply(in_ply, &face_props[0]);
			face_other = get_other_properties_ply(in_ply, offsetof(Face_io, other_props));

			/* grab all the face elements */
			for (j = 0; j < elem_count; j++) {
				Face_io face;
				get_element_ply(in_ply, (void *)&face);

				if (face.nverts != 3) {
					fprintf(stderr, "Face has %d vertices (should be three).\n",
						face.nverts);
					exit(-1);
				}

				/* copy info from the "face" structure */
				tlist[j] = new Triangle;
				tlist[j]->nverts = 3;
				tlist[j]->verts[0] = (Vertex *)face.verts[0];
				tlist[j]->verts[1] = (Vertex *)face.verts[1];
				tlist[j]->verts[2] = (Vertex *)face.verts[2];
				tlist[j]->other_props = face.other_props;
			}
		}
		else
			get_other_element_ply(in_ply);
	}

	/* close the file */
	close_ply(in_ply);

	/* fix up vertex pointers in triangles */
	for (i = 0; i < ntris; i++) {
		tlist[i]->verts[0] = vlist[(int)tlist[i]->verts[0]];
		tlist[i]->verts[1] = vlist[(int)tlist[i]->verts[1]];
		tlist[i]->verts[2] = vlist[(int)tlist[i]->verts[2]];
	}

	/* get rid of triangles that use the same vertex more than once */

	for (i = ntris - 1; i >= 0; i--) {

		Triangle *tri = tlist[i];
		Vertex *v0 = tri->verts[0];
		Vertex *v1 = tri->verts[1];
		Vertex *v2 = tri->verts[2];

		if (v0 == v1 || v1 == v2 || v2 == v0) {
			free(tlist[i]);
			ntris--;
			tlist[i] = tlist[ntris];
		}
	}
}


/******************************************************************************
Write out a polyhedron to a file.
******************************************************************************/

void Polyhedron::write_file(FILE *file)
{
	int i;
	PlyFile *ply;
	char **elist;
	int num_elem_types;

	/*** Write out the transformed PLY object ***/

	elist = get_element_list_ply(in_ply, &num_elem_types);
	ply = write_ply(file, num_elem_types, elist, in_ply->file_type);

	/* describe what properties go into the vertex elements */

	describe_element_ply(ply, "vertex", nverts);
	describe_property_ply(ply, &vert_props[0]);
	describe_property_ply(ply, &vert_props[1]);
	describe_property_ply(ply, &vert_props[2]);
	//  describe_other_properties_ply (ply, vert_other, offsetof(Vertex_io,other_props));

	describe_element_ply(ply, "face", ntris);
	describe_property_ply(ply, &face_props[0]);

	//  describe_other_properties_ply (ply, face_other,
	//                                offsetof(Face_io,other_props));

	//  describe_other_elements_ply (ply, in_ply->other_elems);

	copy_comments_ply(ply, in_ply);
	char mm[1024];
	sprintf(mm, "modified by learnply");
	//  append_comment_ply (ply, "modified by simvizply %f");
	append_comment_ply(ply, mm);
	copy_obj_info_ply(ply, in_ply);

	header_complete_ply(ply);

	/* set up and write the vertex elements */
	put_element_setup_ply(ply, "vertex");
	for (i = 0; i < nverts; i++) {
		Vertex_io vert;

		/* copy info to the "vert" structure */
		vert.x = vlist[i]->x;
		vert.y = vlist[i]->y;
		vert.z = vlist[i]->z;
		vert.other_props = vlist[i]->other_props;

		put_element_ply(ply, (void *)&vert);
	}

	/* index all the vertices */
	for (i = 0; i < nverts; i++)
		vlist[i]->index = i;

	/* set up and write the face elements */
	put_element_setup_ply(ply, "face");

	Face_io face;
	face.verts = new int[3];

	for (i = 0; i < ntris; i++) {

		/* copy info to the "face" structure */
		face.nverts = 3;
		face.verts[0] = tlist[i]->verts[0]->index;
		face.verts[1] = tlist[i]->verts[1]->index;
		face.verts[2] = tlist[i]->verts[2]->index;
		face.other_props = tlist[i]->other_props;

		put_element_ply(ply, (void *)&face);
	}
	put_other_elements_ply(ply);

	close_ply(ply);
	free_ply(ply);
}

void Polyhedron::initialize() {
	icVector3 v1, v2;

	create_pointers();
	calc_edge_length();
	seed = -1;
}

void Polyhedron::finalize() {
	int i;

	for (i = 0; i<ntris; i++) {
		free(tlist[i]->other_props);
		free(tlist[i]);
	}
	for (i = 0; i<nedges; i++) {
		free(elist[i]->tris);
		free(elist[i]);
	}
	for (i = 0; i<nverts; i++) {
		free(vlist[i]->tris);
		free(vlist[i]->other_props);
		free(vlist[i]);
	}

	free(tlist);
	free(elist);
	free(vlist);
	if (!vert_other)
		free(vert_other);
	if (!face_other)
		free(face_other);
}

/******************************************************************************
Find out if there is another face that shares an edge with a given face.

Entry:
f1    - face that we're looking to share with
v1,v2 - two vertices of f1 that define edge

Exit:
return the matching face, or NULL if there is no such face
******************************************************************************/

Triangle *Polyhedron::find_common_edge(Triangle *f1, Vertex *v1, Vertex *v2)
{
	int i, j;
	Triangle *f2;
	Triangle *adjacent = NULL;

	/* look through all faces of the first vertex */

	for (i = 0; i < v1->ntris; i++) {
		f2 = v1->tris[i];
		if (f2 == f1)
			continue;
		/* examine the vertices of the face for a match with the second vertex */
		for (j = 0; j < f2->nverts; j++) {

			/* look for a match */
			if (f2->verts[j] == v2) {

#if 0
				/* watch out for triple edges */

				if (adjacent != NULL) {

					fprintf(stderr, "model has triple edges\n");

					fprintf(stderr, "face 1: ");
					for (k = 0; k < f1->nverts; k++)
						fprintf(stderr, "%d ", f1->iverts[k]);
					fprintf(stderr, "\nface 2: ");
					for (k = 0; k < f2->nverts; k++)
						fprintf(stderr, "%d ", f2->iverts[k]);
					fprintf(stderr, "\nface 3: ");
					for (k = 0; k < adjacent->nverts; k++)
						fprintf(stderr, "%d ", adjacent->iverts[k]);
					fprintf(stderr, "\n");

				}

				/* if we've got a match, remember this face */
				adjacent = f2;
#endif

#if 1
				/* if we've got a match, return this face */
				return (f2);
#endif

			}
		}
	}

	return (adjacent);
}


/******************************************************************************
Create an edge.

Entry:
v1,v2 - two vertices of f1 that define edge
******************************************************************************/

void Polyhedron::create_edge(Vertex *v1, Vertex *v2)
{
	int i, j;
	Triangle *f;

	/* make sure there is enough room for a new edge */

	if (nedges >= max_edges) {

		max_edges += 100;
		Edge **list = new Edge *[max_edges];

		/* copy the old list to the new one */
		for (i = 0; i < nedges; i++)
			list[i] = elist[i];

		/* replace list */
		free(elist);
		elist = list;
	}

	/* create the edge */

	elist[nedges] = new Edge;
	Edge *e = elist[nedges];
	e->index = nedges;
	e->verts[0] = v1;
	e->verts[1] = v2;
	nedges++;

	/* count all triangles that will share the edge, and do this */
	/* by looking through all faces of the first vertex */

	e->ntris = 0;

	for (i = 0; i < v1->ntris; i++) {
		f = v1->tris[i];
		/* examine the vertices of the face for a match with the second vertex */
		for (j = 0; j < 3; j++) {
			/* look for a match */
			if (f->verts[j] == v2) {
				e->ntris++;
				break;
			}
		}
	}

	/* make room for the face pointers (at least two) */
	if (e->ntris < 2)
		e->tris = new Triangle *[2];
	else
		e->tris = new Triangle *[e->ntris];

	/* create pointers from edges to faces and vice-versa */

	e->ntris = 0; /* start this out at zero again for creating ptrs to tris */

	for (i = 0; i < v1->ntris; i++) {

		f = v1->tris[i];

		/* examine the vertices of the face for a match with the second vertex */
		for (j = 0; j < 3; j++)
			if (f->verts[j] == v2) {

				e->tris[e->ntris] = f;
				e->ntris++;

				if (f->verts[(j + 1) % 3] == v1)
					f->edges[j] = e;
				else if (f->verts[(j + 2) % 3] == v1)
					f->edges[(j + 2) % 3] = e;
				else {
					fprintf(stderr, "Non-recoverable inconsistancy in create_edge()\n");
					exit(-1);
				}

				break;  /* we'll only find one instance of v2 */
			}

	}
}


/******************************************************************************
Create edges.
******************************************************************************/

void Polyhedron::create_edges()
{
	int i, j;
	Triangle *f;
	Vertex *v1, *v2;
	double count = 0;

	/* count up how many edges we may require */

	for (i = 0; i < ntris; i++) {
		f = tlist[i];
		for (j = 0; j < f->nverts; j++) {
			v1 = f->verts[j];
			v2 = f->verts[(j + 1) % f->nverts];
			Triangle *result = find_common_edge(f, v1, v2);
			if (result)
				count += 0.5;
			else
				count += 1;
		}
	}

	/*
	printf ("counted %f edges\n", count);
	*/

	/* create space for edge list */

	max_edges = (int)(count + 10);  /* leave some room for expansion */
	elist = new Edge *[max_edges];
	nedges = 0;

	/* zero out all the pointers from faces to edges */

	for (i = 0; i < ntris; i++)
		for (j = 0; j < 3; j++)
			tlist[i]->edges[j] = NULL;

	/* create all the edges by examining all the triangles */

	for (i = 0; i < ntris; i++) {
		f = tlist[i];
		for (j = 0; j < 3; j++) {
			/* skip over edges that we've already created */
			if (f->edges[j])
				continue;
			v1 = f->verts[j];
			v2 = f->verts[(j + 1) % f->nverts];
			create_edge(v1, v2);
		}
	}
}


/******************************************************************************
Create pointers from vertices to faces.
******************************************************************************/

void Polyhedron::vertex_to_tri_ptrs()
{
	int i, j;
	Triangle *f;
	Vertex *v;

	/* zero the count of number of pointers to faces */

	for (i = 0; i < nverts; i++)
		vlist[i]->max_tris = 0;

	/* first just count all the face pointers needed for each vertex */

	for (i = 0; i < ntris; i++) {
		f = tlist[i];
		for (j = 0; j < f->nverts; j++)
			f->verts[j]->max_tris++;
	}

	/* allocate memory for face pointers of vertices */

	for (i = 0; i < nverts; i++) {
		vlist[i]->tris = (Triangle **)
			malloc(sizeof(Triangle *) * vlist[i]->max_tris);
		vlist[i]->ntris = 0;
	}

	/* now actually create the face pointers */

	for (i = 0; i < ntris; i++) {
		f = tlist[i];
		for (j = 0; j < f->nverts; j++) {
			v = f->verts[j];
			v->tris[v->ntris] = f;
			v->ntris++;
		}
	}
}


/******************************************************************************
Find the other triangle that is incident on an edge, or NULL if there is
no other.
******************************************************************************/

Triangle *Polyhedron::other_triangle(Edge *edge, Triangle *tri)
{
	/* search for any other triangle */

	for (int i = 0; i < edge->ntris; i++)
		if (edge->tris[i] != tri)
			return (edge->tris[i]);

	/* there is no such other triangle if we get here */
	return (NULL);
}


/******************************************************************************
Order the pointers to faces that are around a given vertex.

Entry:
v - vertex whose face list is to be ordered
******************************************************************************/

void Polyhedron::order_vertex_to_tri_ptrs(Vertex *v)
{
	int i, j;
	Triangle *f;
	Triangle *fnext;
	int nf;
	int vindex;
	int boundary;
	int count;

	nf = v->ntris;
	f = v->tris[0];

	/* go backwards (clockwise) around faces that surround a vertex */
	/* to find out if we reach a boundary */

	boundary = 0;

	for (i = 1; i <= nf; i++) {

		/* find reference to v in f */
		vindex = -1;
		for (j = 0; j < f->nverts; j++)
			if (f->verts[j] == v) {
				vindex = j;
				break;
			}

		/* error check */
		if (vindex == -1) {
			fprintf(stderr, "can't find vertex #1\n");
			exit(-1);
		}

		/* corresponding face is the previous one around v */
		fnext = other_triangle(f->edges[vindex], f);

		/* see if we've reached a boundary, and if so then place the */
		/* current face in the first position of the vertice's face list */

		if (fnext == NULL) {
			/* find reference to f in v */
			for (j = 0; j < v->ntris; j++)
				if (v->tris[j] == f) {
					v->tris[j] = v->tris[0];
					v->tris[0] = f;
					break;
				}
			boundary = 1;
			break;
		}

		f = fnext;
	}

	/* now walk around the faces in the forward direction and place */
	/* them in order */

	f = v->tris[0];
	count = 0;

	for (i = 1; i < nf; i++) {

		/* find reference to vertex in f */
		vindex = -1;
		for (j = 0; j < f->nverts; j++)
			if (f->verts[(j + 1) % f->nverts] == v) {
				vindex = j;
				break;
			}

		/* error check */
		if (vindex == -1) {
			fprintf(stderr, "can't find vertex #2\n");
			exit(-1);
		}

		/* corresponding face is next one around v */
		fnext = other_triangle(f->edges[vindex], f);

		/* break out of loop if we've reached a boundary */
		count = i;
		if (fnext == NULL) {
			break;
		}

		/* swap the next face into its proper place in the face list */
		for (j = 0; j < v->ntris; j++)
			if (v->tris[j] == fnext) {
				v->tris[j] = v->tris[i];
				v->tris[i] = fnext;
				break;
			}

		f = fnext;
	}
}


/******************************************************************************
Find the index to a given vertex in the list of vertices of a given face.

Entry:
f - face whose vertex list is to be searched
v - vertex to return reference to

Exit:
returns index in face's list, or -1 if vertex not found
******************************************************************************/

int Polyhedron::face_to_vertex_ref(Triangle *f, Vertex *v)
{
	int j;
	int vindex = -1;

	for (j = 0; j < f->nverts; j++)
		if (f->verts[j] == v) {
			vindex = j;
			break;
		}

	return (vindex);
}

/******************************************************************************
Create various face and vertex pointers.
******************************************************************************/

void Polyhedron::create_pointers()
{
	int i;

	/* index the vertices and triangles */

	for (i = 0; i < nverts; i++)
		vlist[i]->index = i;

	for (i = 0; i < ntris; i++)
		tlist[i]->index = i;

	/* create pointers from vertices to triangles */
	vertex_to_tri_ptrs();

	/* make edges */
	create_edges();


	/* order the pointers from vertices to faces */
	for (i = 0; i < nverts; i++) {
		//		if (i %1000 == 0)
		//			fprintf(stderr, "ordering %d of %d vertices\n", i, nverts);
		order_vertex_to_tri_ptrs(vlist[i]);

	}
	/* index the edges */

	for (i = 0; i < nedges; i++) {
		//		if (i %1000 == 0)
		//			fprintf(stderr, "indexing %d of %d edges\n", i, nedges);
		elist[i]->index = i;
	}

}

void Polyhedron::calc_bounding_sphere()
{
	unsigned int i;
	icVector3 min, max;

	for (i = 0; i<nverts; i++) {
		if (i == 0) {
			min.set(vlist[i]->x, vlist[i]->y, vlist[i]->z);
			max.set(vlist[i]->x, vlist[i]->y, vlist[i]->z);
		}
		else {
			if (vlist[i]->x < min.entry[0])
				min.entry[0] = vlist[i]->x;
			if (vlist[i]->x > max.entry[0])
				max.entry[0] = vlist[i]->x;
			if (vlist[i]->y < min.entry[1])
				min.entry[1] = vlist[i]->y;
			if (vlist[i]->y > max.entry[1])
				max.entry[1] = vlist[i]->y;
			if (vlist[i]->z < min.entry[2])
				min.entry[2] = vlist[i]->z;
			if (vlist[i]->z > max.entry[2])
				max.entry[2] = vlist[i]->z;
		}
	}
	center = (min + max) * 0.5;
	radius = length(center - min);
}

void Polyhedron::calc_edge_length()
{
	int i;
	icVector3 v1, v2;

	for (i = 0; i<nedges; i++) {
		v1.set(elist[i]->verts[0]->x, elist[i]->verts[0]->y, elist[i]->verts[0]->z);
		v2.set(elist[i]->verts[1]->x, elist[i]->verts[1]->y, elist[i]->verts[1]->z);
		elist[i]->length = length(v1 - v2);
	}
}

void Polyhedron::calc_face_normals_and_area()
{
	unsigned int i, j;
	icVector3 v0, v1, v2;
	Triangle *temp_t;
	double length[3];

	area = 0.0;
	for (i = 0; i<ntris; i++) {
		for (j = 0; j<3; j++)
			length[j] = tlist[i]->edges[j]->length;
		double temp_s = (length[0] + length[1] + length[2]) / 2.0;
		tlist[i]->area = sqrt(temp_s*(temp_s - length[0])*(temp_s - length[1])*(temp_s - length[2]));

		area += tlist[i]->area;
		temp_t = tlist[i];
		v1.set(vlist[tlist[i]->verts[0]->index]->x, vlist[tlist[i]->verts[0]->index]->y, vlist[tlist[i]->verts[0]->index]->z);
		v2.set(vlist[tlist[i]->verts[1]->index]->x, vlist[tlist[i]->verts[1]->index]->y, vlist[tlist[i]->verts[1]->index]->z);
		v0.set(vlist[tlist[i]->verts[2]->index]->x, vlist[tlist[i]->verts[2]->index]->y, vlist[tlist[i]->verts[2]->index]->z);
		tlist[i]->normal = cross(v0 - v1, v2 - v1);
		normalize(tlist[i]->normal);
	}

	double signedvolume = 0.0;
	icVector3 test = center;
	for (i = 0; i<ntris; i++) {
		icVector3 cent(vlist[tlist[i]->verts[0]->index]->x, vlist[tlist[i]->verts[0]->index]->y, vlist[tlist[i]->verts[0]->index]->z);
		signedvolume += dot(test - cent, tlist[i]->normal)*tlist[i]->area;
	}
	signedvolume /= area;
	if (signedvolume<0)
		orientation = 0;
	else {
		orientation = 1;
		for (i = 0; i<ntris; i++)
			tlist[i]->normal *= -1.0;
	}
}

void sort(unsigned int *A, unsigned int *B, unsigned int *C, unsigned int sid, unsigned int eid) {
	unsigned int i;
	unsigned int *tempA, *tempB, *tempC;
	unsigned int current1, current2, current0;

	if (sid >= eid)
		return;
	sort(A, B, C, sid, (sid + eid) / 2);
	sort(A, B, C, (sid + eid) / 2 + 1, eid);
	tempA = (unsigned int *)malloc(sizeof(unsigned int)*(eid - sid + 1));
	tempB = (unsigned int *)malloc(sizeof(unsigned int)*(eid - sid + 1));
	tempC = (unsigned int *)malloc(sizeof(unsigned int)*(eid - sid + 1));
	for (i = 0; i<eid - sid + 1; i++) {
		tempA[i] = A[i + sid];
		tempB[i] = B[i + sid];
		tempC[i] = C[i + sid];
	}
	current1 = sid;
	current2 = (sid + eid) / 2 + 1;
	current0 = sid;
	while ((current1 <= (sid + eid) / 2) && (current2 <= eid)) {
		if (tempA[current1 - sid] < tempA[current2 - sid]) {
			A[current0] = tempA[current1 - sid];
			B[current0] = tempB[current1 - sid];
			C[current0] = tempC[current1 - sid];
			current1++;
		}
		else if (tempA[current1 - sid] > tempA[current2 - sid]) {
			A[current0] = tempA[current2 - sid];
			B[current0] = tempB[current2 - sid];
			C[current0] = tempC[current2 - sid];
			current2++;
		}
		else {
			if (tempB[current1 - sid] < tempB[current2 - sid]) {
				A[current0] = tempA[current1 - sid];
				B[current0] = tempB[current1 - sid];
				C[current0] = tempC[current1 - sid];
				current1++;
			}
			else {
				A[current0] = tempA[current2 - sid];
				B[current0] = tempB[current2 - sid];
				C[current0] = tempC[current2 - sid];
				current2++;
			}
		}
		current0++;
	}
	if (current1 <= (sid + eid) / 2) {
		for (i = current1; i <= (sid + eid) / 2; i++) {
			A[current0] = tempA[i - sid];
			B[current0] = tempB[i - sid];
			C[current0] = tempC[i - sid];
			current0++;
		}
	}
	if (current2 <= eid) {
		for (i = current2; i <= eid; i++) {
			A[current0] = tempA[i - sid];
			B[current0] = tempB[i - sid];
			C[current0] = tempC[i - sid];
			current0++;
		}
	}

	free(tempA);
	free(tempB);
	free(tempC);
}

void init(void) {
	/* select clearing color */

	glClearColor(0.0, 0.0, 0.0, 0.0);  // background
	glShadeModel(GL_FLAT);
	glPolygonMode(GL_FRONT, GL_FILL);
	//glPolygonMode(GL_FRONT, GL_LINE);

	glDisable(GL_DITHER);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	// may need it
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glEnable(GL_NORMALIZE);
	if (orientation == 0)
		glFrontFace(GL_CW);
	else
		glFrontFace(GL_CCW);
}

void keyboard(unsigned char key, int x, int y) {
	int i;

	Lflag = true;

	/* set escape key to exit */
	switch (key) {
	case 27:
		poly->finalize();  // finalize_everything
		exit(0);
		break;

	case '1':
		display_mode = 0;
		operationNum = 4;
		display();
		break;

	case '2':
		display_mode = 0;
		operationNum = 5;
		display();
		break;

	case '3':
		display_mode = 0;
		operationNum = 6;
		display();
		break;

	case 'x':
		switch (ACSIZE) {
		case 1:
			ACSIZE = 16;
			break;

		case 16:
			ACSIZE = 1;
			break;

		default:
			ACSIZE = 1;
			break;
		}
		fprintf(stderr, "ACSIZE=%d\n", ACSIZE);
		display();
		break;

	case '|':
		this_file = fopen("rotmat.txt", "w");
		for (i = 0; i<4; i++)
			fprintf(this_file, "%f %f %f %f\n", rotmat[i][0], rotmat[i][1], rotmat[i][2], rotmat[i][3]);
		fclose(this_file);
		break;

	case '^':
		this_file = fopen("rotmat.txt", "r");
		for (i = 0; i<4; i++)
			fscanf(this_file, "%f %f %f %f ", (&rotmat[i][0]), (&rotmat[i][1]), (&rotmat[i][2]), (&rotmat[i][3]));
		fclose(this_file);
		display();
		break;

	}
}

Polyhedron::Polyhedron()
{
	nverts = nedges = ntris = 0;
	max_verts = max_tris = 50;

	vlist = new Vertex *[max_verts];
	tlist = new Triangle *[max_tris];
}


void multmatrix(const Matrix m)
{
	int i, j, index = 0;

	GLfloat mat[16];

	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			mat[index++] = m[i][j];

	glMultMatrixf(mat);
}

void set_view(GLenum mode, Polyhedron *poly)
{
	icVector3 up, ray, view;
	GLfloat light_ambient0[] = { 0.3, 0.3, 0.3, 1.0 };
	GLfloat light_diffuse0[] = { 0.7, 0.7, 0.7, 1.0 };
	GLfloat light_specular0[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_ambient1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_diffuse1[] = { 0.5, 0.5, 0.5, 1.0 };
	GLfloat light_specular1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_ambient2[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_diffuse2[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_specular2[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_position[] = { 0.0, 0.0, 0.0, 1.0 };

	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient0);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse0);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular0);
	glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient1);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse1);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular1);


	glMatrixMode(GL_PROJECTION);
	if (mode == GL_RENDER)
		glLoadIdentity();

	if (view_mode == 0)
		glOrtho(-radius_factor, radius_factor, -radius_factor, radius_factor, 0.0, 40.0);
	else
		gluPerspective(45.0, 1.0, 0.1, 40.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	light_position[0] = 5.5;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	light_position[0] = -0.1;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT2, GL_POSITION, light_position);
}

void set_scene(GLenum mode, Polyhedron *poly)
{
	glTranslatef(0.0, 0.0, -3.0);
	multmatrix(rotmat);

	glScalef(1.0 / poly->radius, 1.0 / poly->radius, 1.0 / poly->radius);
	glTranslatef(-poly->center.entry[0], -poly->center.entry[1], -poly->center.entry[2]);
}

void motion(int x, int y) {
	float r[4];
	float xsize, ysize, s, t;

	switch (mouse_mode) {
	case -1:

		xsize = (float)win_width;
		ysize = (float)win_height;

		s = (2.0 * x - win_width) / win_width;
		t = (2.0 * (win_height - y) - win_height) / win_height;

		if ((s == s_old) && (t == t_old))
			return;

		mat_to_quat(rotmat, rvec);
		trackball(r, s_old, t_old, s, t);
		add_quats(r, rvec, rvec);
		quat_to_mat(rvec, rotmat);

		s_old = s;
		t_old = t;

		display();
		break;
	}
}

int processHits(GLint hits, GLuint buffer[])
{
	unsigned int i, j;
	GLuint names, *ptr;
	double smallest_depth = 1.0e+20, current_depth;
	int seed_id = -1;
	unsigned char need_to_update;

	printf("hits = %d\n", hits);
	ptr = (GLuint *)buffer;
	for (i = 0; i < hits; i++) {  /* for each hit  */
		need_to_update = 0;
		names = *ptr;
		ptr++;

		current_depth = (double)*ptr / 0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		current_depth = (double)*ptr / 0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		for (j = 0; j < names; j++) {  /* for each name */
			if (need_to_update == 1)
				seed_id = *ptr - 1;
			ptr++;
		}
	}
	printf("triangle id = %d\n", seed_id);
	return seed_id;
}

void mouse(int button, int state, int x, int y) {
	if (button == GLUT_LEFT_BUTTON || button == GLUT_RIGHT_BUTTON) {
		switch (mouse_mode) {
		case -2:  // no action
			if (state == GLUT_DOWN) {
				float xsize = (float)win_width;
				float ysize = (float)win_height;

				float s = (2.0 * x - win_width) / win_width;
				float t = (2.0 * (win_height - y) - win_height) / win_height;

				s_old = s;
				t_old = t;

				mouse_mode = -1;  // down
				mouse_button = button;
				last_x = x;
				last_y = y;
			}
			break;

		default:
			if (state == GLUT_UP) {
				button = -1;
				mouse_mode = -2;
			}
			break;
		}
	}
	else if (button == GLUT_MIDDLE_BUTTON) {
		if (state == GLUT_DOWN) {  // build up the selection feedback mode

			GLuint selectBuf[win_width];
			GLint hits;
			GLint viewport[4];

			glGetIntegerv(GL_VIEWPORT, viewport);

			glSelectBuffer(win_width, selectBuf);
			(void)glRenderMode(GL_SELECT);

			glInitNames();
			glPushName(0);

			glMatrixMode(GL_PROJECTION);
			glPushMatrix();
			glLoadIdentity();
			/*  create 5x5 pixel picking region near cursor location */
			gluPickMatrix((GLdouble)x, (GLdouble)(viewport[3] - y),
				1.0, 1.0, viewport);

			set_view(GL_SELECT, poly);
			glPushMatrix();
			set_scene(GL_SELECT, poly);
			display_shape(GL_SELECT, poly);
			glPopMatrix();
			glFlush();

			hits = glRenderMode(GL_RENDER);
			poly->seed = processHits(hits, selectBuf);
			display();
		}
	}
}

void display_object()
{
	unsigned int i, j;
	Polyhedron *the_patch = poly;
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	//glPolygonMode(GL_FRONT, GL_LINE);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	for (i = 0; i<poly->ntris; i++) {
		Triangle *temp_t = poly->tlist[i];
		glBegin(GL_POLYGON);
		GLfloat mat_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };

		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);

		glColor3f(1.0, 1.0, 1.0);
		glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);
		for (j = 0; j<3; j++) {
			Vertex *temp_v = temp_t->verts[j];
			glVertex3d(temp_v->x, temp_v->y, temp_v->z);
		}
		glEnd();
	}
}

void display_shape(GLenum mode, Polyhedron *this_poly)
{
	unsigned int i, j;
	GLfloat mat_diffuse[4];

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);

	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	//glPolygonMode(GL_FRONT, GL_LINE);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);

	for (i = 0; i<this_poly->ntris; i++) {
		if (mode == GL_SELECT)
			glLoadName(i + 1);

		Triangle *temp_t = this_poly->tlist[i];

		switch (display_mode) {
		case 0:
			if (i == this_poly->seed) {
				mat_diffuse[0] = 0.0;
				mat_diffuse[1] = 0.0;
				mat_diffuse[2] = 1.0;
				mat_diffuse[3] = 1.0;
			}
			else {
				mat_diffuse[0] = 1.0;
				mat_diffuse[1] = 1.0;
				mat_diffuse[2] = 0.0;
				mat_diffuse[3] = 1.0;
			}
			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			glBegin(GL_POLYGON);
			for (j = 0; j<3; j++) {

				Vertex *temp_v = temp_t->verts[j];
				glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
				if (i == this_poly->seed)
					glColor3f(0.0, 0.0, 1.0);
				else
					glColor3f(1.0, 1.0, 0.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			//generatePolyIDMap(this_poly);
			break;

		case 6:
			glBegin(GL_POLYGON);
			for (j = 0; j<3; j++) {
				Vertex *temp_v = temp_t->verts[j];
				glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);
				glColor3f(1.0, 1.0, 1.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;

		case 10:
			glBegin(GL_POLYGON);
			for (j = 0; j<3; j++) {
				mat_diffuse[0] = 1.0;
				mat_diffuse[1] = 0.0;
				mat_diffuse[2] = 0.0;
				mat_diffuse[3] = 1.0;

				glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);

				Vertex *temp_v = temp_t->verts[j];
				glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);

				glColor3f(1.0, 0.0, 0.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;
		}
	}
	glDisable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
}

void display(void)
{
	GLint viewport[4];
	int jitter;

	glClearColor(1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting
	glGetIntegerv(GL_VIEWPORT, viewport);

	glClear(GL_ACCUM_BUFFER_BIT);
	for (jitter = 0; jitter < ACSIZE; jitter++) {
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		set_view(GL_RENDER, poly);
		glPushMatrix();;
		switch (ACSIZE) {
		case 1:
			glTranslatef(ji1[jitter].x*2.0 / viewport[2], ji1[jitter].y*2.0 / viewport[3], 0.0);
			break;

		case 16:
			glTranslatef(ji16[jitter].x*2.0 / viewport[2], ji16[jitter].y*2.0 / viewport[3], 0.0);
			break;

		default:
			glTranslatef(ji1[jitter].x*2.0 / viewport[2], ji1[jitter].y*2.0 / viewport[3], 0.0);
			break;
		}
		set_scene(GL_RENDER, poly);

		switch (operationNum)
		{
		case 4:
			drawBoundingBox(poly);
			display_shape(GL_RENDER, poly);
			break;

		case 5:
			drawBoundingBox(poly);
			display_shape(GL_RENDER, poly);
			break;

		case 6:
			drawBoundingBox(poly);
			display_shape(GL_RENDER, poly);
			break;

		case -1:
			display_shape(GL_RENDER, poly);
			break;
		}



		glPopMatrix();
		glAccum(GL_ACCUM, 1.0 / ACSIZE);
	}
	glAccum(GL_RETURN, 1.0);
	glFlush();
	glutSwapBuffers();
	glFinish();
}

void Polyhedron::average_normals()
{
	int i, j;

	for (i = 0; i<nverts; i++) {
		vlist[i]->normal = icVector3(0.0);
		for (j = 0; j<vlist[i]->ntris; j++)
			vlist[i]->normal += vlist[i]->tris[j]->normal;
		normalize(vlist[i]->normal);
	}
}

void computeMassCenter(Polyhedron *this_poly, icVector3 &center)
{
	icVector3 total;
	int nverts = this_poly->nverts;
	Vertex **vlist = this_poly->vlist;
	for (unsigned int i = 0; i<nverts; i++) {
		if (i == 0)
			total.set(vlist[i]->x, vlist[i]->y, vlist[i]->z);
		else {
			total.entry[0] += vlist[i]->x;
			total.entry[1] += vlist[i]->y;
			total.entry[2] += vlist[i]->z;
		}
	}
	total /= ((double)nverts);
	center = total;
}

void eigVec(icMatrix3x3 &varience, icVector3 &majEigv, icVector3 &medEigv, icVector3 &minEigv)
{
	float **Nmatrix, **eigmatrix, *eigenvec, mineig, medeig, maxeig;
	int nrot, minind, medind, maxind;
	double x, y, z;

	/* allocate matrices/vector */
	Nmatrix = matrix(1, 3, 1, 3);
	eigmatrix = matrix(1, 3, 1, 3);
	eigenvec = vector(1, 3);

	/* fill up the elements of the matrix */
	for (int i = 0;i<3;++i)
		for (int j = 0;j<3;++j)
			Nmatrix[i + 1][j + 1] = varience.entry[i][j];
	jacobi(Nmatrix, 3, eigenvec, eigmatrix, &nrot);

	/* find smallest eigenvalue */
	mineig = eigenvec[1];
	minind = 1;
	for (int i = 2; i <= 3; ++i) {
		if (eigenvec[i] < mineig) {
			mineig = eigenvec[i];
			minind = i;
		}
	}
	maxeig = eigenvec[1];
	maxind = 1;
	for (int i = 2; i <= 3; ++i) {
		if (eigenvec[i] >= maxeig) {
			maxeig = eigenvec[i];
			maxind = i;
		}
	}
	for (int i = 1;i<4;++i)
	{
		if (i == maxind || i == minind)
			continue;
		medeig = eigenvec[i];
		medind = i;
	}

	majEigv.set(eigmatrix[maxind][1], eigmatrix[maxind][2], eigmatrix[maxind][3]);
	normalize(majEigv);
	medEigv.set(eigmatrix[medind][1], eigmatrix[medind][2], eigmatrix[medind][3]);
	normalize(medEigv);
	minEigv.set(eigmatrix[minind][1], eigmatrix[minind][2], eigmatrix[minind][3]);
	normalize(minEigv);
}


void axisBoundingBox(Polyhedron *this_poly, icVector3 &center, icVector3 &majEigv, icVector3 &medEigv, icVector3 &minEigv)
{
	//icVector3 center;
	majEigv.set(1, 0, 0);
	medEigv.set(0, 1, 0);
	minEigv.set(0, 0, 1);
}

void momentBoundingBox(Polyhedron *this_poly, icVector3 &center, icVector3 &majEigv, icVector3 &medEigv, icVector3 &minEigv)
{
	//icVector3 center;
	computeMassCenter(this_poly, center);

	icMatrix3x3 variance(0.0);
	int nverts = this_poly->nverts;
	Vertex **vlist = this_poly->vlist;
	icMatrix3x3 temp;
	double x;
	double y;
	double z;
	for (unsigned int i = 0; i<nverts; i++) {
		x = vlist[i]->x;
		y = vlist[i]->y;
		z = vlist[i]->z;
		temp.set(x*x, x*y, x*z, y*x, y*y, y*z, z*x, z*y, z*z);
		variance += temp;
	}
	eigVec(variance, majEigv, medEigv, minEigv);

}

void normBoundingBox(Polyhedron *this_poly, icVector3 &center, icVector3 &majEigv, icVector3 &medEigv, icVector3 &minEigv)
{
	computeMassCenter(this_poly, center);

	icMatrix3x3 variance(0.0);
	int ntris = this_poly->ntris;
	Triangle **tlist = this_poly->tlist;

	int nverts = this_poly->nverts;
	Vertex **vlist = this_poly->vlist;


	icMatrix3x3 temp(0.0);
	double x;
	double y;
	double z;
	for (unsigned int i = 0; i<nverts; i++) {

		x = vlist[i]->normal.x;
		y = vlist[i]->normal.y;
		z = vlist[i]->normal.z;

		temp.set(x*x, x*y, x*z, y*x, y*y, y*z, z*x, z*y, z*z);
		variance += temp;
	}

	eigVec(variance, majEigv, medEigv, minEigv);
}

void drawLine(icVector3 &x, icVector3 &y, float width, float r, float g, float b)
{
	glDisable(GL_TEXTURE_2D);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glColor3f(r, g, b);

	glLineWidth(width);
	glBegin(GL_LINES);
	glVertex3d(x.x, x.y, x.z);
	glVertex3d(y.x, y.y, y.z);
	glEnd();

	glDisable(GL_BLEND);
	glEnable(GL_TEXTURE_2D);
}


void drawBoundingBox(Polyhedron *this_poly)
{
	icVector3 center, majEigv, medEigv, minEigv;
	if (operationNum == 4)
		axisBoundingBox(this_poly, center, majEigv, medEigv, minEigv);
	if (operationNum == 5)
		momentBoundingBox(this_poly, center, majEigv, medEigv, minEigv);
	if (operationNum == 6)
		normBoundingBox(this_poly, center, majEigv, medEigv, minEigv);
	int nverts = this_poly->nverts;
	Vertex **vlist = this_poly->vlist;
	unsigned int majMaxInd = 0;
	unsigned int majMinInd = 0;
	unsigned int medMaxInd = 0;
	unsigned int medMinInd = 0;
	unsigned int minMaxInd = 0;
	unsigned int minMinInd = 0;
	icVector3 tempVec(vlist[0]->x, vlist[0]->y, vlist[0]->z);
	tempVec -= center;
	double majDotMax = dot(majEigv, tempVec);
	double majDotMin = dot(majEigv, tempVec);
	double medDotMax = dot(medEigv, tempVec);
	double medDotMin = dot(medEigv, tempVec);
	double minDotMax = dot(minEigv, tempVec);
	double minDotMin = dot(minEigv, tempVec);
	double tempDot;
	for (unsigned int i = 1; i<nverts; i++) {
		tempVec.set(vlist[i]->x, vlist[i]->y, vlist[i]->z);
		tempVec -= center;
		tempDot = dot(majEigv, tempVec);
		if (tempDot>majDotMax)
		{
			majDotMax = tempDot;
			majMaxInd = i;
		}
		else if (tempDot<majDotMin)
		{
			majDotMin = tempDot;
			majMinInd = i;
		}

		tempDot = dot(medEigv, tempVec);
		if (tempDot>medDotMax)
		{
			medDotMax = tempDot;
			medMaxInd = i;
		}
		else if (tempDot<medDotMin)
		{
			medDotMin = tempDot;
			medMinInd = i;
		}

		tempDot = dot(minEigv, tempVec);
		if (tempDot>minDotMax)
		{
			minDotMax = tempDot;
			minMaxInd = i;
		}
		else if (tempDot<minDotMin)
		{
			minDotMin = tempDot;
			minMinInd = i;
		}
	}
	icVector3 majEigvMax = majEigv;
	majEigvMax *= majDotMax;
	majEigvMax += center;
	icVector3 majEigvMin = majEigv;
	majEigvMin *= majDotMin;
	majEigvMin += center;

	icVector3 medEigvMax = medEigv;
	medEigvMax *= medDotMax;
	medEigvMax += center;
	icVector3 medEigvMin = medEigv;
	medEigvMin *= medDotMin;
	medEigvMin += center;

	icVector3 minEigvMax = minEigv;
	minEigvMax *= minDotMax;
	minEigvMax += center;
	icVector3 minEigvMin = minEigv;
	minEigvMin *= minDotMin;
	minEigvMin += center;

	icVector3 tempMajMaxToCenter = majEigvMax - center;
	icVector3 tempMajMinToCenter = majEigvMin - center;
	icVector3 tempMedMaxToCenter = medEigvMax - center;
	icVector3 tempMedMinToCenter = medEigvMin - center;
	icVector3 tempMinMaxToCenter = minEigvMax - center;
	icVector3 tempMinMinToCenter = minEigvMin - center;

	icVector3 point1 = majEigvMax + tempMedMaxToCenter;
	point1 += tempMinMaxToCenter;
	icVector3 point2 = majEigvMin + tempMedMaxToCenter;
	point2 += tempMinMaxToCenter;
	drawLine(point1, point2, 1, 1.0, 0.0, 0.0);

	point1 = majEigvMax + tempMedMaxToCenter;
	point1 += tempMinMinToCenter;
	point2 = majEigvMin + tempMedMaxToCenter;
	point2 += tempMinMinToCenter;
	drawLine(point1, point2, 1, 1.0, 0.0, 0.0);

	point1 = majEigvMax + tempMedMinToCenter;
	point1 += tempMinMaxToCenter;
	point2 = majEigvMin + tempMedMinToCenter;
	point2 += tempMinMaxToCenter;
	drawLine(point1, point2, 1, 1.0, 0.0, 0.0);

	point1 = majEigvMax + tempMedMinToCenter;
	point1 += tempMinMinToCenter;
	point2 = majEigvMin + tempMedMinToCenter;
	point2 += tempMinMinToCenter;
	drawLine(point1, point2, 1, 1.0, 0.0, 0.0);

	point1 = medEigvMax + tempMajMaxToCenter;
	point1 += tempMinMaxToCenter;
	point2 = medEigvMin + tempMajMaxToCenter;
	point2 += tempMinMaxToCenter;
	drawLine(point1, point2, 1, 0.0, 1.0, 0.0);

	point1 = medEigvMax + tempMajMaxToCenter;
	point1 += tempMinMinToCenter;
	point2 = medEigvMin + tempMajMaxToCenter;
	point2 += tempMinMinToCenter;
	drawLine(point1, point2, 1, 0.0, 1.0, 0.0);

	point1 = medEigvMax + tempMajMinToCenter;
	point1 += tempMinMaxToCenter;
	point2 = medEigvMin + tempMajMinToCenter;
	point2 += tempMinMaxToCenter;
	drawLine(point1, point2, 1, 0.0, 1.0, 0.0);

	point1 = medEigvMax + tempMajMinToCenter;
	point1 += tempMinMinToCenter;
	point2 = medEigvMin + tempMajMinToCenter;
	point2 += tempMinMinToCenter;
	drawLine(point1, point2, 1, 0.0, 1.0, 0.0);

	point1 = minEigvMax + tempMajMaxToCenter;
	point1 += tempMedMaxToCenter;
	point2 = minEigvMin + tempMajMaxToCenter;
	point2 += tempMedMaxToCenter;
	drawLine(point1, point2, 1, 0.0, 0.0, 1.0);

	point1 = minEigvMax + tempMajMaxToCenter;
	point1 += tempMedMinToCenter;
	point2 = minEigvMin + tempMajMaxToCenter;
	point2 += tempMedMinToCenter;
	drawLine(point1, point2, 1, 0.0, 0.0, 1.0);

	point1 = minEigvMax + tempMajMinToCenter;
	point1 += tempMedMaxToCenter;
	point2 = minEigvMin + tempMajMinToCenter;
	point2 += tempMedMaxToCenter;
	drawLine(point1, point2, 1, 0.0, 0.0, 1.0);

	point1 = minEigvMax + tempMajMinToCenter;
	point1 += tempMedMinToCenter;
	point2 = minEigvMin + tempMajMinToCenter;
	point2 += tempMedMinToCenter;
	drawLine(point1, point2, 1, 0.0, 0.0, 1.0);
}




