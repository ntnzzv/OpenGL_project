/*  Graphics HW-3

	Natan Zaizev
	נתן זייצב

	Eyal Romano
	איל רומנו

*/

//  
// Created by Ran Dror on April 2018
// Copyright (c) Ran Dror. All right reserved.
// Code can not be copied and/or distributed without the express permission of author
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include <stdarg.h>
#include "GL/glut.h" 
#include "glm.h"
#include "MatrixAlgebraLib.h"

//////////////////////////////////////////////////////////////////////
// Grphics Pipeline section
//////////////////////////////////////////////////////////////////////

#define WIN_SIZE 500
#define TEXTURE_SIZE 512
#define CAMERA_DISTANCE_FROM_AXIS_CENTER 10

typedef struct {
	GLfloat point3D[4];
	GLfloat normal[4];
	GLfloat point3DeyeCoordinates[4];
	GLfloat NormalEyeCoordinates[4];
	GLfloat pointScreen[4];
	GLfloat PixelValue;
	GLfloat TextureCoordinates[3];
} Vertex;

enum ProjectionTypeEnum { ORTHOGRAPHIC = 1, PERSPECTIVE };
enum DisplayTypeEnum { FACE_VERTEXES = 11, FACE_COLOR, LIGHTING_FLAT, LIGHTING_GOURARD, LIGHTING_PHONG, TEXTURE, TEXTURE_LIGHTING_PHONG };
enum DisplayNormalEnum { DISPLAY_NORMAL_YES = 21, DISPLAY_NORMAL_NO };

typedef struct {
	GLfloat ModelMinVec[3]; //(left, bottom, near) of a model.
	GLfloat ModelMaxVec[3]; //(right, top, far) of a model.
	GLfloat CameraPos[3];
	GLfloat ModelScale;
	GLfloat ModelTranslateVector[3];
	enum ProjectionTypeEnum ProjectionType;
	enum DisplayTypeEnum DisplayType;
	enum DisplayNormalEnum DisplayNormals;
	GLfloat Lighting_Diffuse;
	GLfloat Lighting_Specular;
	GLfloat Lighting_Ambient;
	GLfloat Lighting_sHininess;
	GLfloat LightPosition[3];
} GuiParamsForYou;

GuiParamsForYou GlobalGuiParamsForYou;

//written for you
void setPixel(GLint x, GLint y, GLfloat r, GLfloat g, GLfloat b);

//you should write
void ModelProcessing();
void VertexProcessing(Vertex* v);
void FaceProcessing(Vertex* v1, Vertex* v2, Vertex* v3, GLfloat FaceColor[3]);
GLfloat LightingEquation(GLfloat point[3], GLfloat PointNormal[3], GLfloat LightPos[3], GLfloat Kd, GLfloat Ks, GLfloat Ka, GLfloat n);
void DrawLineBresenham(GLint x1, GLint y1, GLint x2, GLint y2, GLfloat r, GLfloat g, GLfloat b);
void BerryCentric(Vertex* v1, Vertex* v2, Vertex* v3, GLfloat faceColor[3]);


GLMmodel* model_ptr;
void ClearColorBuffer();
void DisplayColorBuffer();

void GraphicsPipeline()
{
	static GLuint i;
	static GLMgroup* group;
	static GLMtriangle* triangle;
	Vertex v1, v2, v3;
	GLfloat FaceColor[3];

	//calling ModelProcessing every time refreshing screen
	ModelProcessing();

	//call VertexProcessing for every vertrx
	//and then call FaceProcessing for every face
	group = model_ptr->groups;
	srand(0);
	while (group) {
		for (i = 0; i < group->numtriangles; i++) {
			triangle = &(model_ptr->triangles[group->triangles[i]]);

			MatrixCopy(v1.point3D, &model_ptr->vertices[3 * triangle->vindices[0]], 3);
			v1.point3D[3] = 1;
			MatrixCopy(v1.normal, &model_ptr->normals[3 * triangle->nindices[0]], 3);
			v1.normal[3] = 1;
			if (model_ptr->numtexcoords != 0)
				MatrixCopy(v1.TextureCoordinates, &model_ptr->texcoords[2 * triangle->tindices[0]], 2);
			else {
				v1.TextureCoordinates[0] = -1; v1.TextureCoordinates[1] = -1; v1.TextureCoordinates[2] = -1;
			}
			VertexProcessing(&v1);

			MatrixCopy(v2.point3D, &model_ptr->vertices[3 * triangle->vindices[1]], 3);
			v2.point3D[3] = 1;
			MatrixCopy(v2.normal, &model_ptr->normals[3 * triangle->nindices[1]], 3);
			v2.normal[3] = 1;
			if (model_ptr->numtexcoords != 0)
				MatrixCopy(v2.TextureCoordinates, &model_ptr->texcoords[2 * triangle->tindices[1]], 2);
			else {
				v2.TextureCoordinates[0] = -1; v2.TextureCoordinates[1] = -1; v2.TextureCoordinates[2] = -1;
			}
			VertexProcessing(&v2);

			MatrixCopy(v3.point3D, &model_ptr->vertices[3 * triangle->vindices[2]], 3);
			v3.point3D[3] = 1;
			MatrixCopy(v3.normal, &model_ptr->normals[3 * triangle->nindices[2]], 3);
			v3.normal[3] = 1;
			if (model_ptr->numtexcoords != 0)
				MatrixCopy(v3.TextureCoordinates, &model_ptr->texcoords[2 * triangle->tindices[2]], 2);
			else {
				v3.TextureCoordinates[0] = -1; v3.TextureCoordinates[1] = -1; v3.TextureCoordinates[2] = -1;
			}
			VertexProcessing(&v3);

			FaceColor[0] = (GLfloat)rand() / ((GLfloat)RAND_MAX + 1);
			FaceColor[1] = (GLfloat)rand() / ((GLfloat)RAND_MAX + 1);
			FaceColor[2] = (GLfloat)rand() / ((GLfloat)RAND_MAX + 1);
			FaceProcessing(&v1, &v2, &v3, FaceColor);
		}
		group = group->next;
	}

	DisplayColorBuffer();
}

GLfloat Zbuffer[WIN_SIZE][WIN_SIZE];
GLubyte TextureImage[TEXTURE_SIZE][TEXTURE_SIZE][3];
GLfloat Mmodeling[16];
GLfloat Mlookat[16];
GLfloat Mprojection[16];
GLfloat Mviewport[16];

void ModelProcessing()
{
	GLfloat Mscale[16], Mtranslate[16], Maligner[16];
	GLfloat u[3], v[3], w[3];
	GLfloat center[3] = { 0.0, 0.0, 0.0 }, up[3] = { 0.0, 1.0, 0.0 };
	GLfloat left = -1.0, right = 1.0, top = 1.0, bottom = -1.0;
	GLfloat near = CAMERA_DISTANCE_FROM_AXIS_CENTER - 1.0, far = CAMERA_DISTANCE_FROM_AXIS_CENTER + 1.0;
	GLfloat cx = 0.0, cy = 0.0;



	// ex2-3: calculating translate transformation matrix
	//////////////////////////////////////////////////////////////////////////////////
	M4x4identity(Mscale);
	M4x4identity(Mtranslate);
	M4x4identity(Mmodeling);

	Mscale[0 + 4 * 0] = GlobalGuiParamsForYou.ModelScale;
	Mscale[1 + 4 * 1] = GlobalGuiParamsForYou.ModelScale;
	Mscale[2 + 4 * 2] = GlobalGuiParamsForYou.ModelScale;


	// ex2-3: calculating scale transformation matrix
	//////////////////////////////////////////////////////////////////////////////////
	Mtranslate[0 + 4 * 3] = GlobalGuiParamsForYou.ModelTranslateVector[0];
	Mtranslate[1 + 4 * 3] = GlobalGuiParamsForYou.ModelTranslateVector[2];
	Mtranslate[2 + 4 * 3] = GlobalGuiParamsForYou.ModelTranslateVector[1];

	M4multiplyM4(Mmodeling, Mtranslate, Mscale);

	// ex2-4: calculating lookat transformation matrix
	//////////////////////////////////////////////////////////////////////////////////
	M4x4identity(Mlookat);
	M4x4identity(Mtranslate);
	M4x4identity(Maligner);

	Mtranslate[0 + 4 * 3] = -GlobalGuiParamsForYou.CameraPos[0];
	Mtranslate[1 + 4 * 3] = -GlobalGuiParamsForYou.CameraPos[1];
	Mtranslate[2 + 4 * 3] = -GlobalGuiParamsForYou.CameraPos[2];

	Vminus(w, GlobalGuiParamsForYou.CameraPos, center, 3);
	V3Normalize(w);
	V3cross(u, up, w);
	V3Normalize(u);
	V3cross(v, w, u);

	Maligner[0 + 4 * 0] = u[0];	Maligner[0 + 4 * 1] = u[1];	Maligner[0 + 4 * 2] = u[2];

	Maligner[1 + 4 * 0] = v[0];	Maligner[1 + 4 * 1] = v[1];	Maligner[1 + 4 * 2] = v[2];

	Maligner[2 + 4 * 0] = w[0];	Maligner[2 + 4 * 1] = w[1];	Maligner[2 + 4 * 2] = w[2];

	M4multiplyM4(Mlookat, Maligner, Mtranslate);



	// ex2-2: calculating Orthographic or Perspective projection transformation matrix
	//////////////////////////////////////////////////////////////////////////////////
	M4x4identity(Mprojection);
	M4x4identity(Mscale);
	M4x4identity(Mtranslate);

	if (GlobalGuiParamsForYou.ProjectionType == ORTHOGRAPHIC) {
		
		Mscale[0 + 4 * 0] = (2.0 / (right - left));
		Mscale[1 + 4 * 1] = (2.0 / (top - bottom));
		Mscale[2 + 4 * 2] = (-2.0 / (far - near));

		Mtranslate[0 + 4 * 3] = -1.0 * ((right+left) / (right-left));
		Mtranslate[1 + 4 * 3] = -1.0 * ((top + bottom) / (top - bottom));
		Mtranslate[2 + 4 * 3] = -1.0 * ((far + near) / (far - near));


		M4multiplyM4(Mprojection, Mtranslate, Mscale);
	}

	else {
		Mprojection[0 + 4 * 0] = 2 * near / (right - left);
		Mprojection[1 + 4 * 1] = 2 * near / (top - bottom);
		Mprojection[0 + 4 * 2] = (right + left) / (right - left);
		Mprojection[1 + 4 * 2] = (top + bottom) / (top - bottom);
		Mprojection[2 + 4 * 2] = -((far + near) / (far - near));
		Mprojection[3 + 4 * 2] = -1;
		Mprojection[2 + 4 * 3] = -(2 * far * near) / (far - near);


	}
	// ex2-2: calculating viewport transformation matrix
	//////////////////////////////////////////////////////////////////////////////////
	M4x4identity(Mviewport);
	M4x4identity(Mscale);
	M4x4identity(Mtranslate);

	Mscale[0 + 4 * 0] = cx + WIN_SIZE / 2.0;
	Mscale[1 + 4 * 1] = cy + WIN_SIZE / 2.0;
	Mscale[2 + 4 * 2] = 0.5;

	Mtranslate[0 + 4 * 3] = WIN_SIZE / 2.0;
	Mtranslate[1 + 4 * 3] = WIN_SIZE / 2.0;
	Mtranslate[2 + 4 * 3] = 0.5;

	M4multiplyM4(Mviewport, Mtranslate, Mscale);


	// ex3: clearing color and Z-buffer
	//////////////////////////////////////////////////////////////////////////////////
	ClearColorBuffer(); // setting color buffer to background color
	for (int x = 0; x < WIN_SIZE; x++)
	{
		for (int y = 0; y < WIN_SIZE; y++)
		{
			Zbuffer[x][y] = 1;
		}
	}

}


void VertexProcessing(Vertex* v)
{
	GLfloat point3DafterModelingTrans[4];
	GLfloat temp1[4], temp2[4];
	GLfloat point3D_plusNormal_screen[4];
	GLfloat Mmodeling3x3[9], Mlookat3x3[9];

	// ex2-3: modeling transformation v->point3D --> point3DafterModelingTrans
	//////////////////////////////////////////////////////////////////////////////////
	M4multiplyV4(point3DafterModelingTrans, Mmodeling, v->point3D);



	// ex2-4: lookat transformation point3DafterModelingTrans --> v->point3DeyeCoordinates
	//////////////////////////////////////////////////////////////////////////////////
	//MatrixCopy(v->point3DeyeCoordinates, point3DafterModelingTrans, 4);
M4multiplyV4(v->point3DeyeCoordinates, Mlookat, point3DafterModelingTrans);

// ex2-2: transformation from eye coordinates to screen coordinates v->point3DeyeCoordinates --> v->pointScreen
//////////////////////////////////////////////////////////////////////////////////
//MatrixCopy(v->pointScreen, v->point3DeyeCoordinates, 4);
M4multiplyV4(temp1, Mprojection, v->point3DeyeCoordinates);
M4multiplyV4(v->pointScreen, Mviewport, temp1);


// ex2-5: transformation normal from object coordinates to eye coordinates v->normal --> v->NormalEyeCoordinates
//////////////////////////////////////////////////////////////////////////////////
M3fromM4(Mmodeling3x3, Mmodeling);
M3fromM4(Mlookat3x3, Mlookat);
M3multiplyV3(temp1, Mmodeling3x3, v->normal);
M3multiplyV3(v->NormalEyeCoordinates, Mlookat3x3, temp1);
V3Normalize(v->NormalEyeCoordinates);
v->NormalEyeCoordinates[3] = 1;

// ex2-5: drawing normals 
//////////////////////////////////////////////////////////////////////////////////
if (GlobalGuiParamsForYou.DisplayNormals == DISPLAY_NORMAL_YES) {
	V4HomogeneousDivide(v->point3DeyeCoordinates);
	VscalarMultiply(temp1, v->NormalEyeCoordinates, 0.05, 3);
	Vplus(temp2, v->point3DeyeCoordinates, temp1, 4);
	temp2[3] = 1;
	M4multiplyV4(temp1, Mprojection, temp2);
	M4multiplyV4(point3D_plusNormal_screen, Mviewport, temp1);
	V4HomogeneousDivide(point3D_plusNormal_screen);
	V4HomogeneousDivide(v->pointScreen);
	DrawLineBresenham(round(v->pointScreen[0]), round(v->pointScreen[1]), round(point3D_plusNormal_screen[0]), round(point3D_plusNormal_screen[1]), 0, 0, 1);
}

// ex3: calculating lighting for vertex
//////////////////////////////////////////////////////////////////////////////////

	v->PixelValue = LightingEquation(
		v->point3DeyeCoordinates,
		v->NormalEyeCoordinates,
		GlobalGuiParamsForYou.LightPosition,
		GlobalGuiParamsForYou.Lighting_Diffuse,
		GlobalGuiParamsForYou.Lighting_Specular,
		GlobalGuiParamsForYou.Lighting_Ambient,
		GlobalGuiParamsForYou.Lighting_sHininess);

}



void FaceProcessing(Vertex* v1, Vertex* v2, Vertex* v3, GLfloat FaceColor[3])
{

	V4HomogeneousDivide(v1->pointScreen);
	V4HomogeneousDivide(v2->pointScreen);
	V4HomogeneousDivide(v3->pointScreen);

	if (GlobalGuiParamsForYou.DisplayType == FACE_VERTEXES)
	{
		DrawLineBresenham(round(v1->pointScreen[0]), round(v1->pointScreen[1]), round(v2->pointScreen[0]), round(v2->pointScreen[1]), 1, 1, 1);
		DrawLineBresenham(round(v2->pointScreen[0]), round(v2->pointScreen[1]), round(v3->pointScreen[0]), round(v3->pointScreen[1]), 1, 1, 1);
		DrawLineBresenham(round(v3->pointScreen[0]), round(v3->pointScreen[1]), round(v1->pointScreen[0]), round(v1->pointScreen[1]), 1, 1, 1);
	}
	else {
		//ex3: Barycentric Coordinates and lighting
		//////////////////////////////////////////////////////////////////////////////////

		GLfloat x_min = min(min(v1->pointScreen[0], v2->pointScreen[0]), v3->pointScreen[0]);
		GLfloat x_max = max(max(v1->pointScreen[0], v2->pointScreen[0]), v3->pointScreen[0]);
		GLfloat y_min = min(min(v1->pointScreen[1], v2->pointScreen[1]), v3->pointScreen[1]);
		GLfloat y_max = max(max(v1->pointScreen[1], v2->pointScreen[1]), v3->pointScreen[1]);
		GLfloat pixelDepth, flat_shading_color, gourard_shading_color, phong_shading_color,phong_position[3],phong_normal[3];
		GLfloat A1, B1, C1, A2, B2, C2, A3, B3, C3, x1, y1, x2, y2, x3, y3, alpha, beta, gamma;
		GLint x, y;

		x1 = round(v1->pointScreen[0]);		y1 = round(v1->pointScreen[1]);
		x2 = round(v2->pointScreen[0]);		y2 = round(v2->pointScreen[1]);
		x3 = round(v3->pointScreen[0]);		y3 = round(v3->pointScreen[1]);

		A1 = y1 - y2;						A2 = y2 - y3;						A3 = y3 - y1;
		B1 = x2 - x1;						B2 = x3 - x2;						B3 = x1 - x3;
		C1 = (y2 * x1) - (y1 * x2);			C2 = (y3 * x2) - (y2 * x3);			C3 = (y1 * x3) - (y3 * x1);

	
		for (x = x_min; x < x_max && x >= 0; x++) {
			for (y = y_min; y < y_max&& y >= 0; y++) {

				alpha = (A1 * x + B1 * y + C1) / (A1 * x3 + B1 * y3 + C1);
				beta = (A2 * x + B2 * y + C2) / (A2 * x1 + B2 * y1 + C2);
				gamma = (A3 * x + B3 * y + C3) / (A3 * x2 + B3 * y2 + C3);

				if ((alpha >= 0 && alpha <= 1) && (beta >= 0 && beta <= 1) && (gamma >= 0 && gamma <= 1)) 
				{
					pixelDepth = (alpha * v3->pointScreen[2]) + (beta * v1->pointScreen[2]) + (gamma * v2->pointScreen[2]);

					if (pixelDepth < Zbuffer[x][y])
					{
						Zbuffer[x][y] = pixelDepth;
						
						if (GlobalGuiParamsForYou.DisplayType == FACE_COLOR)
						{
								setPixel(x, y, FaceColor[0], FaceColor[1], FaceColor[2]);
						}
						else if (GlobalGuiParamsForYou.DisplayType == LIGHTING_FLAT)
						{

								flat_shading_color = (v1->PixelValue + v2->PixelValue + v3->PixelValue) / 3.0;
								setPixel(x, y, flat_shading_color, flat_shading_color, flat_shading_color);

						}
						else if (GlobalGuiParamsForYou.DisplayType == LIGHTING_GOURARD)
						{
								gourard_shading_color = (alpha * v3->PixelValue) + (beta * v1->PixelValue) + (gamma * v2->PixelValue);
								setPixel(x, y, gourard_shading_color, gourard_shading_color, gourard_shading_color);
						}
						else if (GlobalGuiParamsForYou.DisplayType == LIGHTING_PHONG)
						{
							for (int i = 0; i < 3; i++)
							{
								phong_position[i] = (alpha * v3->point3DeyeCoordinates[i]) + (beta * v1->point3DeyeCoordinates[i]) + (gamma * v2->point3DeyeCoordinates[i]);
								phong_normal[i] = (alpha * v3->NormalEyeCoordinates[i]) + (beta * v1->NormalEyeCoordinates[i]) + (gamma * v2->NormalEyeCoordinates[i]);
							}
							phong_shading_color = LightingEquation(phong_position, phong_normal, GlobalGuiParamsForYou.LightPosition, GlobalGuiParamsForYou.Lighting_Diffuse, GlobalGuiParamsForYou.Lighting_Specular, GlobalGuiParamsForYou.Lighting_Ambient, GlobalGuiParamsForYou.Lighting_sHininess);
							setPixel(x, y, phong_shading_color, phong_shading_color, phong_shading_color);
						}
						else if (GlobalGuiParamsForYou.DisplayType == TEXTURE) {

							GLfloat s, t, Texture_R, Texture_G, Texture_B, adj_t, adj_s;

							s = (alpha * v3->TextureCoordinates[0]) + (beta * v1->TextureCoordinates[0]) + (gamma * v2->TextureCoordinates[0]);
							t = (alpha * v3->TextureCoordinates[1]) + (beta * v1->TextureCoordinates[1]) + (gamma * v2->TextureCoordinates[1]);
							
							adj_t = round(t * TEXTURE_SIZE);
							adj_s = round(s * TEXTURE_SIZE);

							Texture_R = TextureImage[(int)adj_t][(int)adj_s][0];
							Texture_G = TextureImage[(int)adj_t][(int)adj_s][1];
							Texture_B = TextureImage[(int)adj_t][(int)adj_s][2];
							
							setPixel(x, y,Texture_R / 255,Texture_G / 255,Texture_B / 255);
						}
						else if (GlobalGuiParamsForYou.DisplayType == TEXTURE_LIGHTING_PHONG) {
							GLfloat s, t, Texture_R, Texture_G, Texture_B, adj_t, adj_s;

							s = (alpha * v3->TextureCoordinates[0]) + (beta * v1->TextureCoordinates[0]) + (gamma * v2->TextureCoordinates[0]);
							t = (alpha * v3->TextureCoordinates[1]) + (beta * v1->TextureCoordinates[1]) + (gamma * v2->TextureCoordinates[1]);

							adj_t = round(t * TEXTURE_SIZE);
							adj_s = round(s * TEXTURE_SIZE);

							Texture_R = TextureImage[(int)adj_t][(int)adj_s][0];
							Texture_G = TextureImage[(int)adj_t][(int)adj_s][1];
							Texture_B = TextureImage[(int)adj_t][(int)adj_s][2];

							for (int i = 0; i < 3; i++)
							{
								phong_position[i] = (alpha * v3->point3DeyeCoordinates[i]) + (beta * v1->point3DeyeCoordinates[i]) + (gamma * v2->point3DeyeCoordinates[i]);
								phong_normal[i] = (alpha * v3->NormalEyeCoordinates[i]) + (beta * v1->NormalEyeCoordinates[i]) + (gamma * v2->NormalEyeCoordinates[i]);
							}

							phong_shading_color = LightingEquation(phong_position, phong_normal, GlobalGuiParamsForYou.LightPosition, GlobalGuiParamsForYou.Lighting_Diffuse, GlobalGuiParamsForYou.Lighting_Specular, GlobalGuiParamsForYou.Lighting_Ambient, GlobalGuiParamsForYou.Lighting_sHininess);

							setPixel(x, y, (Texture_R / 255)*phong_shading_color, (Texture_G / 255)*phong_shading_color, (Texture_B / 255)*phong_shading_color);

						}
					}
				}
			}
		}
	}
}


void DrawLineBresenham(GLint x1, GLint y1, GLint x2, GLint y2, GLfloat r, GLfloat g, GLfloat b)
{
	float dx, dy, x, y, a, x1_, y1_, x2_, y2_;

	if ((y2 - y1) > -(x2 - x1)) {
		x1_ = x1;
		y1_ = y1;
		x2_ = x2;
		y2_ = y2;
	}
	else
	{
		x1_ = x2;
		y1_ = y2;
		x2_ = x1;
		y2_ = y1;
	}

	dx = x2_ - x1_;
	dy = y2_ - y1_;
	if (fabs(dx) > fabs(dy)) {
		a = dy / dx;
		y = y1_;
		for (x = x1_; x < x2_; x++) {
			setPixel(x, round(y), 1, 1, 1);
			y = y + a;
		}
	}
	else {
		a = dx / dy;
		x = x1_;
		for (y = y1_; y < y2_; y++) {
			setPixel(round(x), y, 1, 1, 1);
			x = x + a;
		}
	}
}




GLfloat LightingEquation(GLfloat point[3], GLfloat PointNormal[3], GLfloat LightPos[3], GLfloat diffuse, GLfloat specular, GLfloat ambient, GLfloat shine)
{
	//ex3: calculate lighting equation

	//////////////////////////////////////////////////////////////////////////////////
	GLfloat V[3], refl_ray[3], Normal[3], L[3];
	GLfloat I_light = 1.0;

	for (int i = 0; i < 3; ++i)
	{
		Normal[i] = PointNormal[i];
		V[i] = PointNormal[i] - point[i];
		L[i] = LightPos[i] - point[i];
	}

	V3Normalize(Normal);
	V3Normalize(V);
	V3Normalize(L);

	for (int i = 0; i < 3; i++)
	{
		if (V3dot(L,Normal) <= 0) {
			refl_ray[i] = 0;
		}
		else {
			refl_ray[i] = (L[i] + (-2.0 * Normal[i] * (V3dot(Normal, L))));
		}
	}
	V3Normalize(refl_ray);

	return I_light * (max(diffuse * V3dot(Normal, L), 0) + max(specular * pow(V3dot(V, refl_ray), shine), 0) + max(ambient, 0));
	

}





//////////////////////////////////////////////////////////////////////
// GUI section
//////////////////////////////////////////////////////////////////////

//function declerations
void InitGuiGlobalParams();
void drawingCB(void);
void reshapeCB(int width, int height);
void keyboardCB(unsigned char key, int x, int y);
void keyboardSpecialCB(int key, int x, int y);
void MouseClickCB(int button, int state, int x, int y);
void MouseMotionCB(int x, int y);
void menuCB(int value);
void TerminationErrorFunc(char* ErrorString);
void LoadModelFile();
void DisplayColorBuffer();
void drawstr(char* FontName, int FontSize, GLuint x, GLuint y, char* format, ...);
void TerminationErrorFunc(char* ErrorString);
GLubyte* readBMP(char* imagepath, int* width, int* height);

enum FileNumberEnum { TEAPOT = 100, TEDDY, PUMPKIN, COW, SIMPLE_PYRAMID, FIRST_EXAMPLE, SIMPLE_3D_EXAMPLE, SPHERE, TRIANGLE, Z_BUFFER_EXAMPLE, TEXTURE_BOX, TEXTURE_TRIANGLE, TEXTURE_BARREL, TEXTURE_SHEEP };
typedef struct {
	enum FileNumberEnum FileNum;
	GLfloat CameraRaduis;
	GLint   CameraAnleHorizontal;
	GLint   CameraAnleVertical;
	GLint   MouseLastPos[2];
} GuiCalculations;

GuiCalculations GlobalGuiCalculations;

GLuint ColorBuffer[WIN_SIZE][WIN_SIZE][3];

int main(int argc, char** argv)
{
	GLint submenu1_id, submenu2_id, submenu3_id, submenu4_id;

	//initizlizing GLUT
	glutInit(&argc, argv);

	//initizlizing GUI globals
	InitGuiGlobalParams();

	//initializing window
	glutInitWindowSize(WIN_SIZE, WIN_SIZE);
	glutInitWindowPosition(900, 100);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	glutCreateWindow("Computer Graphics HW");

	//registering callbacks
	glutDisplayFunc(drawingCB);
	glutReshapeFunc(reshapeCB);
	glutKeyboardFunc(keyboardCB);
	glutSpecialFunc(keyboardSpecialCB);
	glutMouseFunc(MouseClickCB);
	glutMotionFunc(MouseMotionCB);

	//registering and creating menu
	submenu1_id = glutCreateMenu(menuCB);
	glutAddMenuEntry("open teapot.obj", TEAPOT);
	glutAddMenuEntry("open teddy.obj", TEDDY);
	glutAddMenuEntry("open pumpkin.obj", PUMPKIN);
	glutAddMenuEntry("open cow.obj", COW);
	glutAddMenuEntry("open Simple3Dexample.obj", SIMPLE_3D_EXAMPLE);
	glutAddMenuEntry("open SimplePyramid.obj", SIMPLE_PYRAMID);
	glutAddMenuEntry("open sphere.obj", SPHERE);
	glutAddMenuEntry("open triangle.obj", TRIANGLE);
	glutAddMenuEntry("open FirstExample.obj", FIRST_EXAMPLE);
	glutAddMenuEntry("open ZbufferExample.obj", Z_BUFFER_EXAMPLE);
	glutAddMenuEntry("open TriangleTexture.obj", TEXTURE_TRIANGLE);
	glutAddMenuEntry("open box.obj", TEXTURE_BOX);
	glutAddMenuEntry("open barrel.obj", TEXTURE_BARREL);
	glutAddMenuEntry("open sheep.obj", TEXTURE_SHEEP);
	submenu2_id = glutCreateMenu(menuCB);
	glutAddMenuEntry("Orthographic", ORTHOGRAPHIC);
	glutAddMenuEntry("Perspective", PERSPECTIVE);
	submenu3_id = glutCreateMenu(menuCB);
	glutAddMenuEntry("Face Vertexes", FACE_VERTEXES);
	glutAddMenuEntry("Face Color", FACE_COLOR);
	glutAddMenuEntry("Lighting Flat", LIGHTING_FLAT);
	glutAddMenuEntry("Lighting Gourard", LIGHTING_GOURARD);
	glutAddMenuEntry("Lighting Phong", LIGHTING_PHONG);
	glutAddMenuEntry("Texture", TEXTURE);
	glutAddMenuEntry("Texture lighting Phong", TEXTURE_LIGHTING_PHONG);
	submenu4_id = glutCreateMenu(menuCB);
	glutAddMenuEntry("Yes", DISPLAY_NORMAL_YES);
	glutAddMenuEntry("No", DISPLAY_NORMAL_NO);
	glutCreateMenu(menuCB);
	glutAddSubMenu("Open Model File", submenu1_id);
	glutAddSubMenu("Projection Type", submenu2_id);
	glutAddSubMenu("Display type", submenu3_id);
	glutAddSubMenu("Display Normals", submenu4_id);
	glutAttachMenu(GLUT_RIGHT_BUTTON);

	LoadModelFile();

	//starting main loop
	glutMainLoop();
}

void drawingCB(void)
{
	GLenum er;

	char DisplayString1[200], DisplayString2[200];

	//clearing the background
	glClearColor(0, 0, 0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//initializing modelview transformation matrix
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	GraphicsPipeline();

	glColor3f(0, 1, 0);
	sprintf(DisplayString1, "Scale:%.1f , Translate: (%.1f,%.1f,%.1f), Camera angles:(%d,%d) position:(%.1f,%.1f,%.1f) ", GlobalGuiParamsForYou.ModelScale, GlobalGuiParamsForYou.ModelTranslateVector[0], GlobalGuiParamsForYou.ModelTranslateVector[1], GlobalGuiParamsForYou.ModelTranslateVector[2], GlobalGuiCalculations.CameraAnleHorizontal, GlobalGuiCalculations.CameraAnleVertical, GlobalGuiParamsForYou.CameraPos[0], GlobalGuiParamsForYou.CameraPos[1], GlobalGuiParamsForYou.CameraPos[2]);
	drawstr("helvetica", 12, 15, 25, DisplayString1);
	sprintf(DisplayString2, "Lighting reflection - Diffuse:%1.2f, Specular:%1.2f, Ambient:%1.2f, sHininess:%1.2f", GlobalGuiParamsForYou.Lighting_Diffuse, GlobalGuiParamsForYou.Lighting_Specular, GlobalGuiParamsForYou.Lighting_Ambient, GlobalGuiParamsForYou.Lighting_sHininess);
	drawstr("helvetica", 12, 15, 10, DisplayString2);

	//swapping buffers and displaying
	glutSwapBuffers();

	//check for errors
	er = glGetError();  //get errors. 0 for no error, find the error codes in: https://www.opengl.org/wiki/OpenGL_Error
	if (er) printf("error: %d\n", er);
}


void LoadModelFile()
{
	int width, height;
	GLubyte* ImageData;

	if (model_ptr) {
		glmDelete(model_ptr);
		model_ptr = 0;
	}

	switch (GlobalGuiCalculations.FileNum) {
	case TEAPOT:
		model_ptr = glmReadOBJ("teapot.obj");
		break;
	case TEDDY:
		model_ptr = glmReadOBJ("teddy.obj");
		break;
	case PUMPKIN:
		model_ptr = glmReadOBJ("pumpkin.obj");
		break;
	case COW:
		model_ptr = glmReadOBJ("cow.obj");
		break;
	case SIMPLE_PYRAMID:
		model_ptr = glmReadOBJ("SimplePyramid.obj");
		break;
	case FIRST_EXAMPLE:
		model_ptr = glmReadOBJ("FirstExample.obj");
		break;
	case SIMPLE_3D_EXAMPLE:
		model_ptr = glmReadOBJ("Simple3Dexample.obj");
		break;
	case SPHERE:
		model_ptr = glmReadOBJ("sphere.obj");
		break;
	case TRIANGLE:
		model_ptr = glmReadOBJ("triangle.obj");
		break;
	case Z_BUFFER_EXAMPLE:
		model_ptr = glmReadOBJ("ZbufferExample.obj");
		break;
	case TEXTURE_TRIANGLE:
		model_ptr = glmReadOBJ("TriangleTexture.obj");
		ImageData = readBMP("TriangleTexture.bmp", &width, &height);
		if (width != TEXTURE_SIZE || height != TEXTURE_SIZE)
			TerminationErrorFunc("Invalid texture size");
		memcpy(TextureImage, ImageData, TEXTURE_SIZE * TEXTURE_SIZE * 3);
		free(ImageData);
		break;
	case TEXTURE_BOX:
		model_ptr = glmReadOBJ("box.obj");
		ImageData = readBMP("box.bmp", &width, &height);
		if (width != TEXTURE_SIZE || height != TEXTURE_SIZE)
			TerminationErrorFunc("Invalid texture size");
		memcpy(TextureImage, ImageData, TEXTURE_SIZE * TEXTURE_SIZE * 3);
		free(ImageData);
		break;
	case TEXTURE_BARREL:
		model_ptr = glmReadOBJ("barrel.obj");
		ImageData = readBMP("barrel.bmp", &width, &height);
		if (width != TEXTURE_SIZE || height != TEXTURE_SIZE)
			TerminationErrorFunc("Invalid texture size");
		memcpy(TextureImage, ImageData, TEXTURE_SIZE * TEXTURE_SIZE * 3);
		free(ImageData);
		break;
	case TEXTURE_SHEEP:
		model_ptr = glmReadOBJ("sheep.obj");
		ImageData = readBMP("sheep.bmp", &width, &height);
		if (width != TEXTURE_SIZE || height != TEXTURE_SIZE)
			TerminationErrorFunc("Invalid texture size");
		memcpy(TextureImage, ImageData, TEXTURE_SIZE * TEXTURE_SIZE * 3);
		free(ImageData);
		break;
	default:
		TerminationErrorFunc("File number not valid");
		break;
	}

	if (!model_ptr)
		TerminationErrorFunc("can not load 3D model");
	//glmUnitize(model_ptr);  //"unitize" a model by translating it

	//to the origin and scaling it to fit in a unit cube around
	//the origin
	glmFacetNormals(model_ptr);  //adding facet normals
	glmVertexNormals(model_ptr, 90.0);  //adding vertex normals

	glmBoundingBox(model_ptr, GlobalGuiParamsForYou.ModelMinVec, GlobalGuiParamsForYou.ModelMaxVec);
}

void ClearColorBuffer()
{
	GLuint x, y;
	for (y = 0; y < WIN_SIZE; y++) {
		for (x = 0; x < WIN_SIZE; x++) {
			ColorBuffer[y][x][0] = 0;
			ColorBuffer[y][x][1] = 0;
			ColorBuffer[y][x][2] = 0;
		}
	}
}

void setPixel(GLint x, GLint y, GLfloat r, GLfloat g, GLfloat b)
{
	if (x >= 0 && x < WIN_SIZE && y >= 0 && y < WIN_SIZE) {
		ColorBuffer[y][x][0] = round(r * 255);
		ColorBuffer[y][x][1] = round(g * 255);
		ColorBuffer[y][x][2] = round(b * 255);
	}
}

void DisplayColorBuffer()
{
	GLuint x, y;
	glBegin(GL_POINTS);
	for (y = 0; y < WIN_SIZE; y++) {
		for (x = 0; x < WIN_SIZE; x++) {
			glColor3ub(min(255, ColorBuffer[y][x][0]), min(255, ColorBuffer[y][x][1]), min(255, ColorBuffer[y][x][2]));
			glVertex2f(x + 0.5, y + 0.5);   // The 0.5 is to target pixel
		}
	}
	glEnd();
}


void InitGuiGlobalParams()
{
	GlobalGuiCalculations.FileNum = TEAPOT;
	GlobalGuiCalculations.CameraRaduis = CAMERA_DISTANCE_FROM_AXIS_CENTER;
	GlobalGuiCalculations.CameraAnleHorizontal = 0;
	GlobalGuiCalculations.CameraAnleVertical = 0;

	GlobalGuiParamsForYou.CameraPos[0] = 0;
	GlobalGuiParamsForYou.CameraPos[1] = 0;
	GlobalGuiParamsForYou.CameraPos[2] = GlobalGuiCalculations.CameraRaduis;

	GlobalGuiParamsForYou.ModelScale = 1;

	GlobalGuiParamsForYou.ModelTranslateVector[0] = 0;
	GlobalGuiParamsForYou.ModelTranslateVector[1] = 0;
	GlobalGuiParamsForYou.ModelTranslateVector[2] = 0;
	GlobalGuiParamsForYou.DisplayType = FACE_VERTEXES;
	GlobalGuiParamsForYou.ProjectionType = ORTHOGRAPHIC;
	GlobalGuiParamsForYou.DisplayNormals = DISPLAY_NORMAL_NO;
	GlobalGuiParamsForYou.Lighting_Diffuse = 0.75;
	GlobalGuiParamsForYou.Lighting_Specular = 0.2;
	GlobalGuiParamsForYou.Lighting_Ambient = 0.2;
	GlobalGuiParamsForYou.Lighting_sHininess = 40;
	GlobalGuiParamsForYou.LightPosition[0] = 10;
	GlobalGuiParamsForYou.LightPosition[1] = 5;
	GlobalGuiParamsForYou.LightPosition[2] = 0;

}


void reshapeCB(int width, int height)
{
	if (width != WIN_SIZE || height != WIN_SIZE)
	{
		glutReshapeWindow(WIN_SIZE, WIN_SIZE);
	}

	//update viewport
	glViewport(0, 0, width, height);

	//clear the transformation matrices (load identity)
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	//projection
	gluOrtho2D(0, WIN_SIZE, 0, WIN_SIZE);
}


void keyboardCB(unsigned char key, int x, int y) {
	switch (key) {
	case 27:
		exit(0);
		break;
	case '+':
		GlobalGuiParamsForYou.ModelScale += 0.1;
		glutPostRedisplay();
		break;
	case '-':
		GlobalGuiParamsForYou.ModelScale -= 0.1;
		glutPostRedisplay();
		break;
	case 'd':
	case 'D':
		GlobalGuiParamsForYou.Lighting_Diffuse += 0.05;
		glutPostRedisplay();
		break;
	case 'c':
	case 'C':
		GlobalGuiParamsForYou.Lighting_Diffuse -= 0.05;
		glutPostRedisplay();
		break;
	case 's':
	case 'S':
		GlobalGuiParamsForYou.Lighting_Specular += 0.05;
		glutPostRedisplay();
		break;
	case 'x':
	case 'X':
		GlobalGuiParamsForYou.Lighting_Specular -= 0.05;
		glutPostRedisplay();
		break;
	case 'a':
	case 'A':
		GlobalGuiParamsForYou.Lighting_Ambient += 0.05;
		glutPostRedisplay();
		break;
	case 'z':
	case 'Z':
		GlobalGuiParamsForYou.Lighting_Ambient -= 0.05;
		glutPostRedisplay();
		break;
	case 'h':
	case 'H':
		GlobalGuiParamsForYou.Lighting_sHininess += 1;
		glutPostRedisplay();
		break;
	case 'n':
	case 'N':
		GlobalGuiParamsForYou.Lighting_sHininess -= 1;
		glutPostRedisplay();
		break;
	default:
		printf("Key not valid (language shold be english)\n");
	}
}


void keyboardSpecialCB(int key, int x, int y)
{
	switch (key) {
	case GLUT_KEY_LEFT:
		GlobalGuiParamsForYou.ModelTranslateVector[0] -= 0.1;
		glutPostRedisplay();
		break;
	case GLUT_KEY_RIGHT:
		GlobalGuiParamsForYou.ModelTranslateVector[0] += 0.1;
		glutPostRedisplay();
		break;
	case GLUT_KEY_DOWN:
		GlobalGuiParamsForYou.ModelTranslateVector[2] -= 0.1;
		glutPostRedisplay();
		break;
	case GLUT_KEY_UP:
		GlobalGuiParamsForYou.ModelTranslateVector[2] += 0.1;
		glutPostRedisplay();
		break;
	}
}


void MouseClickCB(int button, int state, int x, int y)
{
	GlobalGuiCalculations.MouseLastPos[0] = x;
	GlobalGuiCalculations.MouseLastPos[1] = y;
}

void MouseMotionCB(int x, int y)
{
	GlobalGuiCalculations.CameraAnleHorizontal += (x - GlobalGuiCalculations.MouseLastPos[0]) / 40;
	GlobalGuiCalculations.CameraAnleVertical -= (y - GlobalGuiCalculations.MouseLastPos[1]) / 40;

	if (GlobalGuiCalculations.CameraAnleVertical > 30)
		GlobalGuiCalculations.CameraAnleVertical = 30;
	if (GlobalGuiCalculations.CameraAnleVertical < -30)
		GlobalGuiCalculations.CameraAnleVertical = -30;

	GlobalGuiCalculations.CameraAnleHorizontal = (GlobalGuiCalculations.CameraAnleHorizontal + 360) % 360;
	//	GlobalGuiCalculations.CameraAnleVertical   = (GlobalGuiCalculations.CameraAnleVertical   + 360) % 360;

	GlobalGuiParamsForYou.CameraPos[0] = GlobalGuiCalculations.CameraRaduis * sin((float)(GlobalGuiCalculations.CameraAnleVertical + 90) * M_PI / 180) * cos((float)(GlobalGuiCalculations.CameraAnleHorizontal + 90) * M_PI / 180);
	GlobalGuiParamsForYou.CameraPos[2] = GlobalGuiCalculations.CameraRaduis * sin((float)(GlobalGuiCalculations.CameraAnleVertical + 90) * M_PI / 180) * sin((float)(GlobalGuiCalculations.CameraAnleHorizontal + 90) * M_PI / 180);
	GlobalGuiParamsForYou.CameraPos[1] = GlobalGuiCalculations.CameraRaduis * cos((float)(GlobalGuiCalculations.CameraAnleVertical + 90) * M_PI / 180);
	glutPostRedisplay();
}

void menuCB(int value)
{
	switch (value) {
	case ORTHOGRAPHIC:
	case PERSPECTIVE:
		GlobalGuiParamsForYou.ProjectionType = value;
		glutPostRedisplay();
		break;
	case FACE_VERTEXES:
	case FACE_COLOR:
	case LIGHTING_FLAT:
	case LIGHTING_GOURARD:
	case LIGHTING_PHONG:
		GlobalGuiParamsForYou.DisplayType = value;
		glutPostRedisplay();
		break;
	case TEXTURE:
		GlobalGuiParamsForYou.DisplayType = value;
		glutPostRedisplay();
		break;
	case TEXTURE_LIGHTING_PHONG:
		GlobalGuiParamsForYou.DisplayType = value;
		glutPostRedisplay();
		break;
	case DISPLAY_NORMAL_YES:
	case DISPLAY_NORMAL_NO:
		GlobalGuiParamsForYou.DisplayNormals = value;
		glutPostRedisplay();
		break;
	default:
		GlobalGuiCalculations.FileNum = value;
		LoadModelFile();
		glutPostRedisplay();
	}
}



void drawstr(char* FontName, int FontSize, GLuint x, GLuint y, char* format, ...)
{
	va_list args;
	char buffer[255], * s;

	GLvoid* font_style = GLUT_BITMAP_TIMES_ROMAN_10;

	font_style = GLUT_BITMAP_HELVETICA_10;
	if (strcmp(FontName, "helvetica") == 0) {
		if (FontSize == 12)
			font_style = GLUT_BITMAP_HELVETICA_12;
		else if (FontSize == 18)
			font_style = GLUT_BITMAP_HELVETICA_18;
	}
	else if (strcmp(FontName, "times roman") == 0) {
		font_style = GLUT_BITMAP_TIMES_ROMAN_10;
		if (FontSize == 24)
			font_style = GLUT_BITMAP_TIMES_ROMAN_24;
	}
	else if (strcmp(FontName, "8x13") == 0) {
		font_style = GLUT_BITMAP_8_BY_13;
	}
	else if (strcmp(FontName, "9x15") == 0) {
		font_style = GLUT_BITMAP_9_BY_15;
	}

	va_start(args, format);
	vsprintf(buffer, format, args);
	va_end(args);

	glRasterPos2i(x, y);
	for (s = buffer; *s; s++)
		glutBitmapCharacter(font_style, *s);
}


void TerminationErrorFunc(char* ErrorString)
{
	char string[256];
	printf(ErrorString);
	fgets(string, 256, stdin);

	exit(0);
}

GLubyte* readBMP(char* imagepath, int* width, int* height)
{
	unsigned char header[54]; // Each BMP file begins by a 54-bytes header
	unsigned int dataPos;     // Position in the file where the actual data begins
	unsigned int imageSize;   // = width*height*3
	unsigned char* data;
	unsigned char tmp;
	int i;

	// Open the file
	FILE* file = fopen(imagepath, "rb");
	if (!file) {
		TerminationErrorFunc("Image could not be opened\n");
	}

	if (fread(header, 1, 54, file) != 54) { // If not 54 bytes read : problem
		TerminationErrorFunc("Not a correct BMP file\n");
	}

	if (header[0] != 'B' || header[1] != 'M') {
		TerminationErrorFunc("Not a correct BMP file\n");
	}

	// Read ints from the byte array
	dataPos = *(int*) & (header[0x0A]);
	imageSize = *(int*) & (header[0x22]);
	*width = *(int*) & (header[0x12]);
	*height = *(int*) & (header[0x16]);

	// Some BMP files are misformatted, guess missing information
	if (imageSize == 0)
		imageSize = *width * *height * 3; // 3 : one byte for each Red, Green and Blue component
	if (dataPos == 0)
		dataPos = 54; // The BMP header is done that way

	// Create a buffer
	data = malloc(imageSize * sizeof(GLubyte));

	// Read the actual data from the file into the buffer
	fread(data, 1, imageSize, file);


	//swap the r and b values to get RGB (bitmap is BGR)
	for (i = 0; i < *width * *height; i++)
	{
		tmp = data[i * 3];
		data[i * 3] = data[i * 3 + 2];
		data[i * 3 + 2] = tmp;
	}


	//Everything is in memory now, the file can be closed
	fclose(file);

	return data;
}