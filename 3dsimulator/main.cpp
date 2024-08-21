//iisph wird gemacht
//theorie: zweite ableitung vom druckfeld
//scalarprodukt
//omega immer 0.5
// konzeptionell riesige matrix mit jakobi gelöst
// divergenz der geschwindigkeit = 0 = kein dichtefehler
// gleichung 2 : links approkimation: rechts SPH  gleichung 2 wird nicht gelößt, da zu rechenaufwendig und unendlich viele lösungen und geschwindigkeitsäderungen auf symetrischen kräften basieren sollen
//gamma erstmal auf 1
//mittlere gleichung gleichung 14 nachtpartikel, partikel und dann andersrum der kernelfunktion
// omega mal (sf-AP) immer kleiner
// durchschnittlicher fehler, absolutbetrag
// omega * sf/aff statt pl mit null initialisieren
//adaptiver cfl an zeitschritt anpassen
//diagonalelement mit tatsächlicher dichte
//maximalgeschwindigkeit mal zeitschritt soll kleiner als partikelgröße sein
#include"imgui.h"
#include"implot.h"
#include"imgui_impl_glfw.h"
#include"imgui_impl_opengl3.h"
#include "Eigen/Dense"
#include <stdio.h>
#include <stdlib.h>
#include<omp.h>
#include <vector>
#include <algorithm>

#include <GL/glew.h>
#include <GLFW/glfw3.h>
GLFWwindow * window;

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/norm.hpp>
using namespace glm;

#include <controls.hpp>
#include <shader.hpp>
#include <texture.hpp>
#include <stdexcept>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <set>
#include <string>
#include <fstream>
#include <Cell.h>
#include <neighboursearch.h>
#include <scenarios.h>
#include <densitycomputation.h>
#include <rigidbodies.h>
#include <helperz.h>
#include <kernelz.h>
#include <simulationframes.h>
#include <massinits.h>
#include <Affz.h>
#include <boundarypres.h>
#include <iisphparticle.h>
#include<part.h>
#include <chrono>
#include <sstream>
#include <unordered_set>


//particleamount___________________________________________________________________________________________________
int var_nx = 36;
int var_ny = 30;
int var_nz = 30;
int var_oneway = 50;
int var_fluidpart = var_nx * var_ny * var_nz;
int var_MaxParticles = var_nx * var_ny * var_nz + (var_oneway * var_oneway * 12);
int var_spezialboundpart = 0;
//constants__________________________________________________________________
float pi_F = glm::pi<float>();
const float h =0.015625;// 0.0039090625;//
const float p0 = 1000;
const float p0p0 = p0*p0;
const float alpha = 1 / (4 * pi_F * h * h * h); //8 / (pi_F*h*h*h);// 
const float alphaTwoD = 5 / (14 * pi_F * h * h);
const float m = (p0 * h * h *h);
float deltaT = 0.001;
float viscosity = 0.7;
float visc_fluid = 0.7;// 0.0000011;
float visc_boundary = 0.3;// 0.0000011;
float k = 10000;
float cfl_max = 0.4;
float deltaTmax = 0.01;
float cfl = 0;
bool cfl_cond = false;
float avgErr = 0.1;
glm::vec3 gamma = glm::vec3(1, 1, 1);
float gammafloat = 0.7f;
float gammapresup = 0.7f;
float gammadens = 1.f;
float gammapres = 0.7f;
float gammabound = 0.7f;
float omega = 0.5;
float gravity = -9.81f;
//hash___________________________________________________________________________________________
int gridResolution = var_oneway * var_oneway * var_oneway;
float cellsize = 2 * h;
int hashsize = var_MaxParticles;
const int p1 = 73856093;
const int p2 = 19349663;
const int p3 = 83492791;
float searchRadius = 2 * h;
//debug_________________________________________________________________________________________________
std::vector<float> densitys;
std::vector<float> densitysnew;
std::vector<float> currenttimes; 
std::vector<float> convergence;
float currentdens;
float neighbouraverage = 0;
float maxvel = 0.f;
float maxdist = 0.f;
//scenarios________________________________________________________________________________________________
bool start = false;
bool tesla = false;
bool teslaclosed = false;
bool realtesla = false;
bool watercol = false;
bool rotatinganimation = false;
bool dambreaktestscen = false;
bool smalldambreaktestscen = false;
bool smallwatercol = false;
bool TwoDwatercol = true;
bool TwoDDambreak = false;

//__________________________________________________________________________________________________
float sizefac = 6;
float distfac = 1;
float overallmaxvel = 0.f;
float cloloroffset = 9.f;
int observing = var_MaxParticles + 1;
float numNeighofobserving;
//importexportstuf_____________________________________________________________________________________________
bool exportdata = false;
bool exportimage = false;
bool exportanimation = false;
bool exportscene = false;
bool importscene = false;
bool exportdens = false;
//__________________________________________________________________________________________________________
glm::vec3 eindrittelvector = glm::vec3(1 / 3, 1 / 3, 1 / 3);
int currentiter = 0;
int currentitermax = 500;
float denistyerrormax =  0.1;
bool newiisph = false;
bool iisph = false;
bool ssph = false;
bool TwoDsimul = true;
int numexport = 0;
float animationtime = 0.f;
float neighSearchTime = 0.f;
float drawingTime = 0.f;
float othercomputationTime = 0.f;
float kerneltime = 0.f;
float densitytime = 0.f;
float nonpresatime = 0.f;
float predveltime = 0.f;
float computeallsftime = 0.f;
float computeallafftime = 0.f;
float looptime = 0.f;
float everythingtime = 0.f;
float overheadtime = 0.f;
bool makesingle = true;
float teslatime = 0.f;
bool maketesla = false;
bool maketeslaclosed = false;
int itertesla = 0;
int boundarya = 10;
int obstaclea = 10;
const int overallmaxpart = 18*(50 * 50 * 50 + (55 * 55 * 12));
float denserrold = 0.f;
float maxavgdensdeviation = 0.f;
bool singlewall = true;
int avgiter= 0;
int maxiter = 0;
float avgcomptime = 0;
float maxcomptime = 0;
int animationstep = 0;
int watercolheight = 50;
float totaltimeanimated = 0;
float simcomptime= 0; 
float max_singledens = 0;
int num_rot_part = 140;
float baseRotationSpeed = 0.25;
float angle = 0;
float minypos = 0;
float maxypos = 0;
int animatedparticles = 0;
int animationneighbourscount = 60;
char changename[128];
std::string firstdatapath;
glm::vec3 oldpos;
bool showallpart = true;
bool showboundandouter = false;
bool showouter = false;
bool absinterrupt = false;
bool fixeddt = false; 
float maxdenserralltime = 0; 
bool paussimul = true;
float oldanimstep = -1;
float totalComputationTime = 0;
float overallmaxcfl = 0;
float surfacetension = 0;
bool highlightunderpop = false;
bool highlightextremefast = false;
int numofp0high = 0;
float upperviualbord = 100; 
float lowervisualbord = -10;
bool colorvel = true;
bool colordens = false;
bool colorpres = false;
bool setboundmass = true;
bool setpartmass = true;
bool addfloating = false;
float gammapart = 1;
float jitterfac = 0.01;
float maxdenserrold = 0;
bool resetvalues = false;
bool updateviewpos = true;
std::vector<float> maxvels;
std::vector<int> usemefordens;
std::vector<float> cfls;
std::vector<float>iterations;
std::vector < float > densdiffes;
float maxdiffdens = 0;
bool highlightusedforcomp = false;
bool highlightdenserr = false;
float gridbreite = 20;
float gridhöhe = 0;
float densobserving = 0;
float massobserving = 0;
float nonpresaobserving = 0;
float presaobserving = 0;
float sumkernelderobserving = 0;
float sumkernelderobservingy = 0;
float sumkernelderobservingz = 0;
float sumkernelobserving = 0;
float allrigidmass = 0;
float observingposx = 0;
float observingposy = 0;
float observingposz = 0;
bool exportconvergence = false;
int currentitermin = 2;
glm::vec3 posofcenterofmass = glm::vec3(0.f, 0.f, 0.f);
glm::vec3 velofcenterofmass = glm::vec3(0.f, 0.f, 0.f);
glm::mat3x3 rotMat = glm::mat3x3(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f);
glm::vec3 xCM = glm::vec3(0.f,0.f,0.f);
glm::vec3 vCM = glm::vec3(0.f, 0.f, 0.f);
glm::mat3 A(0.f);
glm::vec3 L = glm::vec3(0.f, 0.f, 0.f);
glm::mat3 I_inv(0.f);
glm::mat3 inertiaTensor(0.f);
glm::vec3 omegarigidbody = glm::vec3(0.f, 0.f, 0.f);
glm::mat3 inertiaTensorInverse;
glm::vec3 torque(0.f, 0.f, 0.f);
std::vector<std::vector<int>> uniformgidvec1D;
bool uniformgridneighsearch = true;
bool pressuredextrapolatebound = false;
float overallminxpos = 0.f;
float overallminypos = 0.f;
float overallminzpos = 0.f;
float massfac = 1;
float avgsincomptime = 0;
float avgfps = 0;
int simstepcnter = 0;
float fps = 0;
bool exportcfl = false;
bool usewholefluidforcalc = false;
bool usewholefluidtodevide = false;
//_________________________________________________________________________________
int windowedWidth = 2048;//small1600;
int windowedHeight = 1152;//small900;
int* wWidth = &windowedWidth;
int* wHeight = &windowedHeight;
glm::vec3 position = glm::vec3(15, 1, 50);
//initlists__________________________________________________________________
//iisphparticle ParticlesContainer[MaxParticles];
//Particle ParticlesContainer[MaxParticles];
std::unordered_map<int, Cell> hashmap;
std::vector<iisphparticle> var_PartC(var_MaxParticles);

void StepSimulation(std::vector<iisphparticle>& var_PartC, std::unordered_map<int, Cell>& hashmap) {
	if (iisph) {
		var_iisph(var_PartC, hashmap);
	}
	if (ssph) {
		ssphAlgo(var_PartC, hashmap);
	}
	if (TwoDsimul) {
		if (ssph) {
			ssphAlgotwoD(var_PartC, hashmap);
		}
		else {
			twoDiisph(var_PartC, hashmap);
		}
	}
	animationstep++;
	currenttimes.push_back(totaltimeanimated +deltaT);
	if (exportanimation) {
		animationtime += deltaT;
	}
}
int main(void)
{
	// Initialize GLFWc
	if (!glfwInit())
	{
		fprintf(stderr, "Failed to initialize GLFW\n");
		getchar();
		return -1;
	}

	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make macOS happy; should not be needed
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	// Open a window and create its OpenGL context
	window = glfwCreateWindow(windowedWidth, windowedHeight, "crazy great fluid simulation", NULL, NULL);
	if (window == NULL) {
		fprintf(stderr, "Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible. Try the 2.1 version of the tutorials.\n");
		getchar();
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);

	// Initialize GLEW
	glewExperimental = true; // Needed for core profile
	if (glewInit() != GLEW_OK) {
		fprintf(stderr, "Failed to initialize GLEW\n");
		getchar();
		glfwTerminate();
		return -1;
	}
	// Ensure we can capture the escape key being pressed below
	glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
	// Hide the mouse and enable unlimited movement
	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);

	// Set the mouse at the center of the screen
	glfwPollEvents();
	glfwSetCursorPos(window, windowedWidth / 2, windowedHeight / 2);

	// Dark blue background
	glClearColor(1.0f, 1.0f, 1.f, 0.0f);
	// Enable depth test
	glEnable(GL_DEPTH_TEST);
	// Accept fragment if it is closer to the camera than the former one
	glDepthFunc(GL_LESS);

	GLuint VertexArrayID;
	glGenVertexArrays(1, &VertexArrayID);
	glBindVertexArray(VertexArrayID);


	// Create and compile our GLSL program from the shaders
	GLuint programID = LoadShaders("Particle.vertexshader", "Particle.fragmentshader");

	// Vertex shader
	GLuint CameraRight_worldspace_ID = glGetUniformLocation(programID, "CameraRight_worldspace");
	GLuint CameraUp_worldspace_ID = glGetUniformLocation(programID, "CameraUp_worldspace");
	GLuint ViewProjMatrixID = glGetUniformLocation(programID, "VP");

	// fragment shader
	GLuint TextureID = glGetUniformLocation(programID, "myTextureSampler");
	/*
	std::vector<GLfloat> g_particule_position_size_data(50 * 50 * 50 + (55 * 55 * 12));
	std::vector<GLubyte> g_particule_color_data(50 * 50 * 50 + (55 * 55 * 12));
	*/
	static GLfloat* g_particule_position_size_data = new GLfloat[overallmaxpart];
	static GLubyte* g_particule_color_data = new GLubyte[overallmaxpart];
	

	//mystuff
	computeMatricesFromInputs();
	glm::mat4 ProjectionMatrix = getProjectionMatrix();
	glm::mat4 ViewMatrix = getViewMatrix();

	// We will need the camera's position in order to sort the particles
	// w.r.t the camera's distance.
	// There should be a getCameraPosition() function in common/controls.cpp, 
	// but this works too.
	glm::vec3 CameraPosition(glm::inverse(ViewMatrix)[3]);

	glm::mat4 ViewProjectionMatrix = ProjectionMatrix * ViewMatrix;
	GLuint Texture = loadDDS("particle.DDS");

	// The VBO containing the 4 vertices of the particles.
	// Thanks to instancing, they will be shared by all particles.
	static const GLfloat g_vertex_buffer_data[] = {
		 -0.5f, -0.5f, 0.0f,
		  0.5f, -0.5f, 0.0f,
		 -0.5f,  0.5f, 0.0f,
		  0.5f,  0.5f, 0.0f,
	};
	GLuint billboard_vertex_buffer;
	glGenBuffers(1, &billboard_vertex_buffer);
	glBindBuffer(GL_ARRAY_BUFFER, billboard_vertex_buffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(g_vertex_buffer_data), g_vertex_buffer_data, GL_STATIC_DRAW);

	// The VBO containing the positions and sizes of the particles
	GLuint particles_position_buffer;
	glGenBuffers(1, &particles_position_buffer);
	glBindBuffer(GL_ARRAY_BUFFER, particles_position_buffer);
	// Initialize with empty (NULL) buffer : it will be updated later, each frame.
	glBufferData(GL_ARRAY_BUFFER, overallmaxpart * 4 * sizeof(GLfloat), NULL, GL_STREAM_DRAW);

	// The VBO containing the colors of the particles
	GLuint particles_color_buffer;
	glGenBuffers(1, &particles_color_buffer);
	glBindBuffer(GL_ARRAY_BUFFER, particles_color_buffer);
	// Initialize with empty (NULL) buffer : it will be updated later, each frame.
	glBufferData(GL_ARRAY_BUFFER, overallmaxpart * 4 * sizeof(GLubyte), NULL, GL_STREAM_DRAW);
	//__________________________________________________________________________________________________________________________________________________________________________________
	double lastTime = glfwGetTime();
	densitysnew.push_back(0.f);
	densitys.push_back(0.f);
	currenttimes.push_back(0.f);

	//imgui__________________________________________________________________________________________________________________________________________________________________________________
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	//ImPlot::CreateContext();
	ImGuiIO& io = ImGui::GetIO(); (void)io;
	ImGui::StyleColorsDark();
	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL3_Init("#version 330");
	
	do
	{
		auto start_alltime = std::chrono::high_resolution_clock::now();

		// Clear the screen
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		//imgui
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();
		//imgui end
		double currentTime = glfwGetTime();
		double delta = currentTime - lastTime;
		lastTime = currentTime;

		
		computeMatricesFromInputs();
		glm::mat4 ProjectionMatrix = getProjectionMatrix();
		glm::mat4 ViewMatrix = getViewMatrix();

		// We will need the camera's position in order to sort the particles
		// w.r.t the camera's distance.
		// There should be a getCameraPosition() function in common/controls.cpp, 
		// but this works too.
		glm::vec3 CameraPosition(glm::inverse(ViewMatrix)[3]);

		glm::mat4 ViewProjectionMatrix = ProjectionMatrix * ViewMatrix;
		//init
		if (start) {
			gridbreite = var_oneway *2;
			gridhöhe = var_oneway *2;
			cellsize = 2 * h;
			uniformgidvec1D.resize((var_oneway *2) * var_oneway *2 * var_oneway *2);
			upperviualbord = (var_oneway+10)*h;
			lowervisualbord = -1*h;
			if (updateviewpos) {
				position = glm::vec3(4, 2, 30);
			}
			
			firstdatapath = "breakingDam";
			firstdatapath.copy(changename, 127);
			var_init(CameraPosition, var_PartC);
			start = false;
			maketesla = false;
			maketeslaclosed = false;
			sizefac = 1 / h * 0.1;
			if (setboundmass) {
				makeboundmass(var_PartC, hashmap);
			}
			if (setpartmass) {
				makepartmass(var_PartC, hashmap);
			}
		}
		if (dambreaktestscen) {
			gridbreite = 250;
			gridhöhe = 200;
			cellsize = 2 * h;
			uniformgidvec1D.resize((200) * gridbreite * gridhöhe);
			upperviualbord = 58*h;
			lowervisualbord = -1*h;
			fixeddt = false;
			currentitermax = 100;
			denistyerrormax = 0.1;
			visc_fluid = 0.025;
			visc_boundary = 0.0000025;
			absinterrupt = false;
			if (updateviewpos) {
				position = glm::vec3(14, 5, 50);
			}
			
			firstdatapath = "breakingDamTest";
			firstdatapath.copy(changename, 127);
			exportanimation = false;
			gammafloat = 1;
			deltaTmax = 0.001;
			cfl_max = 0.79;
			deltaT = 0.02;
			dambreaktest(CameraPosition, var_PartC);
			dambreaktestscen = false;
			maketesla = false;
			maketeslaclosed = false;
			if (setboundmass) {
				makeboundmass(var_PartC, hashmap);
			}
			if (setpartmass) {
				makepartmass(var_PartC, hashmap);
			}
		}if (smalldambreaktestscen) {
			gridbreite = 250;
			gridhöhe = 200;
			cellsize = 2 * h;
			uniformgidvec1D.resize((200) * gridbreite * gridhöhe);
			upperviualbord = 58*h;
			lowervisualbord = -1*h;
			if (updateviewpos) {
				position = glm::vec3(15, 5, 50);
			}
			
			firstdatapath = "breakingDamTestSmall";
			firstdatapath.copy(changename, 127);
			fixeddt = false;
			//surfacetension = 0.08;
			currentitermax = 1000;
			denistyerrormax = 0.1;
			//jitterfac = 0;
			//gammapres = 0.7;
			visc_fluid = 0.025;
			visc_boundary = 0.0000025;
			absinterrupt = true;
			singlewall = true;
			exportanimation = false;
			//gammafloat = 0.7;
			deltaTmax = 0.0009;
			cfl_max = 0.5;
			deltaT = 0.01;
			//gammapresup = 0.7;
			smalldambreaktest(CameraPosition, var_PartC);
			smalldambreaktestscen = false;
			maketesla = false;
			maketeslaclosed = false;
			TwoDsimul = false;
			resetvalues = true;
			iisph = true;
			if (setboundmass) {
				makeboundmass(var_PartC, hashmap);
			}
			//setpartmass = true;
			if (setpartmass) {
				makepartmass(var_PartC, hashmap);
			}
			makefloatingpartmass(var_PartC, hashmap);
			initrigidbodies(var_PartC, hashmap);
		}
		
		if (realtesla) {
			gridbreite = 250;
			gridhöhe = 200;
			cellsize = 2 * h;
			uniformgidvec1D.resize((200)* gridbreite* gridhöhe);
			upperviualbord = 210;
			lowervisualbord = -10;
			obstaclea = 50;
			if (updateviewpos) {
				position = glm::vec3(12, 6, 40);
			}
			firstdatapath = "Tesla";
			firstdatapath.copy(changename, 127);
			gammafloat = 1.0;;
			deltaTmax = 0.01;
			deltaT = 0.01;
			visc_boundary = 0.3;
			visc_fluid = 1;
			//var_teslavalve(CameraPosition, var_PartC);
			var_real_teslavalve(CameraPosition, var_PartC);
			realtesla = false;
			maketesla = true;
			maketeslaclosed = false;
			if (setboundmass) {
				makeboundmass(var_PartC, hashmap);
			}
			if (setpartmass) {
				makepartmass(var_PartC, hashmap);
			}
		}if (tesla) {
			firstdatapath = "Obstacles";
			var_teslavalve(CameraPosition, var_PartC);
			tesla = false;
			maketesla = true;
			maketeslaclosed = false;
			TwoDsimul = false;
			gridbreite = 250;
			gridhöhe = 250;
			cellsize = 2 * h;
			uniformgidvec1D.resize((200) * gridbreite * gridhöhe);
			upperviualbord = 310 * h;
			lowervisualbord = -10 * h;
			obstaclea = 50;
			if (updateviewpos) {
				position = glm::vec3(12, 6, 40);
			}

			firstdatapath.copy(changename, 127);
			gammafloat = 1.0;
			deltaTmax = 0.001;
			deltaT = 0.01;
			visc_boundary = 0.0003;
			visc_fluid = 0.01;
			if (setboundmass) {
				makeboundmass(var_PartC, hashmap);
			}
			if (setpartmass) {
				makepartmass(var_PartC, hashmap);
			}
		}
		if (teslaclosed) {
			firstdatapath = "ObstaclesClosed";
			teslaclosed = false;
			maketeslaclosed = true;
			maketesla = false;
			iisph = true;
			ssph = false;
			TwoDsimul = false;
			var_teslavalveclosed(CameraPosition, var_PartC);

			gridbreite = 250;
			gridhöhe = 250;
			cellsize = 2 * h;
			uniformgidvec1D.resize((200)* gridbreite* gridhöhe);
			upperviualbord = 310*h;
			lowervisualbord = -10*h;
			obstaclea = 50;
			if (updateviewpos) {
				position = glm::vec3(12, 6, 40);
			}

			firstdatapath.copy(changename, 127);
			gammafloat = 1.0;
			deltaTmax = 0.001;
			deltaT = 0.01;
			visc_boundary = 0.0003;
			visc_fluid = 0.01;
			
			if (setboundmass) {
				makeboundmass(var_PartC, hashmap);
			}
			if (setpartmass) {
				makepartmass(var_PartC, hashmap);
			}
		}
		if (watercol) {
			gridbreite = 200;
			gridhöhe = 200;
			cellsize = 2 * h;
			uniformgidvec1D.resize((200 + watercolheight)* gridbreite* gridhöhe);
			gammapres = 0.7;
			deltaTmax = 0.0008;
			upperviualbord = 200 + watercolheight;
			lowervisualbord = -10;
			k = 500000;
			if (updateviewpos) {
				position = glm::vec3(1, 7, 45);
			}
			
			firstdatapath = "WaterColumn";
			firstdatapath.copy(changename, 127);
			exportanimation = false;
			boundarya = 70;
			absinterrupt = true;
			paussimul = true;
			singlewall = true;
			denistyerrormax = 0.1;
			visc_fluid = 0.025;
			visc_boundary = 0.0025;
			watercolumn(CameraPosition, var_PartC);
			watercol = false;
			TwoDsimul = false;
			iisph = true;
			ssph = false;
			if (setboundmass) {
				makeboundmass(var_PartC, hashmap);
			}
			if (setpartmass) {
				makepartmass(var_PartC, hashmap);
			}
			makefloatingpartmass(var_PartC, hashmap);
			initrigidbodies(var_PartC, hashmap);
		}
		if (smallwatercol) {
			gridbreite = 200;
			gridhöhe = 200;
			cellsize = 2 * h;
			uniformgidvec1D.resize((200 + watercolheight)* gridbreite*gridhöhe);
			//gammapres = 0.5;
			//gammabound = 1.3;
			deltaTmax = 0.001;
			upperviualbord = 200 + watercolheight*h;
			lowervisualbord = -1*h;
			k = 5500000;
			if (updateviewpos) {
				position = glm::vec3(1, 7, 45);
			}
			firstdatapath = "SmallWaterColumn";
			firstdatapath.copy(changename, 127);
			exportanimation = false;
			//boundarya = 70;
			absinterrupt = true;
			paussimul = true;
			//jitterfac = 0;
			gammapres = 0.7;
			//setboundmass = true;
			denistyerrormax = 0.1;
			visc_fluid = 0.025;
			visc_boundary = 0.0000025;
			singlewall = true;
			watercolumnsmall(CameraPosition, var_PartC);
			smallwatercol = false;
			if (setboundmass) {
				makeboundmass(var_PartC, hashmap);
			}
			if (setpartmass) {
				makepartmass(var_PartC, hashmap);
			}
			makefloatingpartmass(var_PartC, hashmap);
			initrigidbodies(var_PartC, hashmap);
			TwoDsimul = false;
			iisph = true;
			ssph = false;
		}
		if (TwoDwatercol) {
			gridbreite = 20;
			gridhöhe = 0;
			cellsize = 2 * h;
			uniformgidvec1D.resize((200+watercolheight)*gridbreite);
			TwoDsimul = true;
			iisph = false; 
			ssph = false; 
			//jitterfac = 0.01f;
			//watercolheight = 10;
			deltaTmax = 0.0009;
			cfl_max = 0.15;
			//jitterfac = 0;
			//gammabound = 0.6;
			//gammapres = 0.5;
			//singlewall = false;
			//setpartmass = true;
			//setboundmass = true;
			upperviualbord = (watercolheight+200)*h;
			lowervisualbord = 0;
			k = 6000000;
			if (updateviewpos) {
				position = glm::vec3(1, 4, 25);
			}
			
			firstdatapath = "/2dWatercol50/gammabound" + std::to_string(gammabound) +"gammapres"+ std::to_string(gammapres);
			firstdatapath.copy(changename, 127);
			exportanimation = false;
			boundarya = 70;
			absinterrupt = true;
			paussimul = true;
			denistyerrormax = 0.1;
			visc_fluid = 0.0025;
			visc_boundary = 0.0000025;
			cloloroffset = 50;
			gammapres = 0.7;
			watercolumnTwoD(CameraPosition, var_PartC);
			
			TwoDwatercol = false;
			if (setboundmass) {
				makeboundmassTwoD(var_PartC, hashmap);
			}
			if (setpartmass) {
				makepartmassTwoD(var_PartC, hashmap);
			}
			makefloatingpartmass2d(var_PartC, hashmap);
			initrigidbodies(var_PartC, hashmap);
			resetvalues = true;
			
		}
		if (TwoDDambreak) {
			TwoDsimul = true;
			iisph = false;
			ssph = false;
			//jitterfac = 0.f;
			gridbreite = 300;
			gridhöhe = 0;
			cellsize = 2 * h;
			uniformgidvec1D.resize((200)* gridbreite);
			watercolheight = 50;
			deltaTmax = 0.0009;
			cfl_max = 0.7;
			var_nx = 20;
			var_ny = 80;
			gammabound = 0.7;
			gammapres = 0.5;
			//setpartmass = true;
			//setboundmass = true;
			upperviualbord =150*h;
			lowervisualbord = -1*h;
			k = 6000000;
			if (updateviewpos) {
				position = glm::vec3(10, 5, 35);
			}
			firstdatapath = "SmallWaterColumn";
			firstdatapath.copy(changename, 127);
			exportanimation = false;
			boundarya = 70;
			absinterrupt = true;
			paussimul = true;
			denistyerrormax = 0.1;
			visc_fluid = 0.0025;
			visc_boundary = 0.0000025;
			cloloroffset = 10;
			DambreaktestTwoD(CameraPosition, var_PartC);
			TwoDDambreak = false;
			if (setboundmass) {
				makeboundmassTwoD(var_PartC, hashmap);
			}
			if (setpartmass) {
				makepartmassTwoD(var_PartC, hashmap);
			}
			makefloatingpartmass2d(var_PartC, hashmap);
			initrigidbodies(var_PartC, hashmap);
		}
		if (rotatinganimation) {
			gridbreite = var_oneway *2;
			gridhöhe = var_oneway *2;
			cellsize = 2 * h;
			uniformgidvec1D.resize((gridbreite)* gridbreite* gridhöhe);
			upperviualbord = 200;
			lowervisualbord = -10;
			if (updateviewpos) {
				position = glm::vec3(1, 2, 30);
			}
			
			firstdatapath = "Blender";
			firstdatapath.copy(changename, 127);
			exportanimation = false;
			num_rot_part = 184;
			var_oneway = 50;
			var_nx = 20;
			var_ny = 49;
			var_nz = 48;
			baseRotationSpeed = 0.0;
			denistyerrormax = 0.1;
			deltaTmax = 0.001;
			upperviualbord = 200 + watercolheight * h;
			lowervisualbord = -10;
			visc_fluid = 0.025;
			visc_boundary = 0.0000025;
			moving_boundary(CameraPosition, var_PartC);
			rotatinganimation = false;
			if (setboundmass) {
				makeboundmass(var_PartC, hashmap);
			}
			if (setpartmass) {
				makepartmass(var_PartC, hashmap);
			}
		}
		if (exportimage) {
			std::string filenamefirst(changename);
			std::string filename("../exportdata/numPart" + std::to_string(var_MaxParticles) + "v" + std::to_string(viscosity) + "k" + std::to_string(k) + filenamefirst +".csv");
			std::fstream file_out;
			file_out.open(filename, std::ios_base::out);
			if (!file_out.is_open()) {
				std::cerr << "Fehler beim Öffnen der Datei." << std::endl;
				return -1;
			}
			for (const auto& part : var_PartC) {
				if (part.isboundary == false && part.pos.y < upperviualbord && part.pos.y > lowervisualbord) {
					std::string xPos = std::to_string(part.pos.x);
					std::string yPos = std::to_string(part.pos.y);
					std::string zPos = std::to_string(part.pos.z);
					std::string normv = std::to_string(glm::length(part.vel));
					file_out << xPos << ";" << yPos << ";" << zPos << ";" << normv << "\n";
				}
			}
			file_out.flush();
			file_out.close();
			std::cout << "finito" << std::endl;
			std::cout << "S key pressed!" << std::endl;
			exportimage = false;
		}
		if (exportscene) {
			std::string filenamefirst(changename);
			std::cout << firstdatapath << std::endl;
			std::string filename("../exportdata/Scene"+ filenamefirst +".csv");
			std::fstream file_out;
			file_out.open(filename, std::ios_base::out);
			if (!file_out.is_open()) {
				std::cerr << "Fehler beim Öffnen der Datei." << std::endl;
				return -1;
			}
			for (const auto& part : var_PartC) {
				std::string xPos = std::to_string(part.pos.x);
				std::string yPos = std::to_string(part.pos.y);
				std::string zPos = std::to_string(part.pos.z);
				std::string xAcc = std::to_string(part.acc.x);
				std::string yAcc = std::to_string(part.acc.y);
				std::string zAcc = std::to_string(part.acc.z);
				std::string xVel = std::to_string(part.vel.x);
				std::string yVel = std::to_string(part.vel.y);
				std::string zVel = std::to_string(part.vel.z);
				std::string bound = std::to_string(part.isboundary);
				file_out << bound << ";" << xPos << ";" << yPos << ";" << zPos << ";" << xAcc << ";" << yAcc << ";" << zAcc << "; " << xVel << ";" << yVel << ";" << zVel << "\n";
			}
			file_out.flush();
			file_out.close();
			std::cout << "finito" << std::endl;
			std::cout << "Scene exported!" << std::endl;
			exportscene = false;
		}
		if (exportdens) {
			std::string filenamefirst(changename);
			std::cout << firstdatapath << std::endl;
			std::string filename("../exportdata/density" + filenamefirst + ".csv");
			std::fstream file_out;
			file_out.open(filename, std::ios_base::out);
			if (!file_out.is_open()) {
				std::cerr << "Fehler beim Öffnen der Datei." << std::endl;
				return -1;
			}
			for (int i = 0; i < densitys.size(); i++) {
				file_out << currenttimes[i] << ";" << densitysnew[i] << ";" << densitys[i] << "\n";
			}
			file_out.flush();
			file_out.close();
			std::cout << "Density exported!" << std::endl;
			exportdens = false;
		}
		if (exportconvergence) {
			std::string filenamefirst(changename);
			std::cout << firstdatapath << std::endl;
			std::string filename("../exportdata/convergence" + filenamefirst + ".csv");
			std::fstream file_out;
			file_out.open(filename, std::ios_base::out);
			if (!file_out.is_open()) {
				std::cerr << "Fehler beim Öffnen der Datei." << std::endl;
				return -1;
			}
			for (int i = 0; i < convergence.size(); i++) {
				file_out << i << ";" << convergence[i] << "\n";
			}
			file_out.flush();
			file_out.close();
			std::cout << "Convergence exported!" << std::endl;
			exportconvergence = false;
		}
		if (exportcfl) {
			std::string filenamefirst(changename);
			std::cout << firstdatapath << std::endl;
			std::string filename("../exportdata/cfl" + filenamefirst + ".csv");
			std::fstream file_out;
			file_out.open(filename, std::ios_base::out);
			if (!file_out.is_open()) {
				std::cerr << "Fehler beim Öffnen der Datei." << std::endl;
				return -1;
			}
			for (int i = 0; i < cfls.size(); i++) {
				file_out << i << ";" << cfls[i] << "\n";
			}
			file_out.flush();
			file_out.close();
			std::cout << "cfl exported!" << std::endl;
			exportcfl = false;
		}
		if (importscene) {
			std::string filenamefirst(changename);
			std::cout << firstdatapath << std::endl;
			std::string filename("../exportdata/Scene" + filenamefirst + ".csv");
			std::ifstream file(filename);
			std::vector<std::vector<std::string>> data;
			std::string line;

			// Check if the file is open
			if (!file.is_open()) {
				std::cerr << "Could not open the file!" << std::endl;
			}
			// Read the file line by line
			while (std::getline(file, line)) {
				std::stringstream lineStream(line);
				std::string cell;
				std::vector<std::string> row;

				// Read each cell separated by a comma
				while (std::getline(lineStream, cell, ';')) {
					row.push_back(cell);
				}
				// Add the row to the data vector
				data.push_back(row);
			}
			var_fluidpart = 0;
			var_MaxParticles = 0;
			for (std::vector<std::string> roww : data) {
				if (roww[0] == "0") {
					var_fluidpart++;
				}
				var_MaxParticles++;
			}
			hashsize = var_MaxParticles;
			var_PartC.resize(0);
			var_PartC.resize(var_MaxParticles);

			for (int i = 0; i < var_MaxParticles; i++) {
				var_PartC[i].resetvalues();
				var_PartC[i].ismovingboundary = false;
				if (data[i][0] == "0") {
					var_PartC[i].isboundary = false;
					var_PartC[i].a = 250;
				}
				else {
					var_PartC[i].isboundary = true;
					var_PartC[i].a = boundarya;
				}
				var_PartC[i].pos.x = std::stof(data[i][1]);
				var_PartC[i].pos.y = std::stof(data[i][2]);
				var_PartC[i].pos.z = std::stof(data[i][3]);
				var_PartC[i].acc.x = std::stof(data[i][4]);
				var_PartC[i].acc.y = std::stof(data[i][5]);
				var_PartC[i].acc.z = std::stof(data[i][6]);
				var_PartC[i].vel.x = std::stof(data[i][7]);
				var_PartC[i].vel.y = std::stof(data[i][8]);
				var_PartC[i].vel.z = std::stof(data[i][9]);
				var_PartC[i].index = i;
			}
			file.close();
			std::cout << "Scene imported!" << std::endl;
			importscene = false;
			exportanimation = false;
		}
		if ((exportanimation && animationtime >=0.04)) {
			std::string filenamefirst(changename);
			std::string filename("../exportdata/Scene" + filenamefirst + std::to_string(animationneighbourscount) + std::to_string(numexport) + ".csv");
			//std::string filename("../exportdata/"+ firstdatapath + std::to_string(animationneighbourscount)+ std::to_string(numexport) + ".csv");
			numexport++;
			std::fstream file_out;
			file_out.open(filename, std::ios_base::out);
			if (!file_out.is_open()) {
				std::cerr << "Fehler beim Öffnen der Datei." << std::endl;
				return -1;
			}
//#pragma omp parallel for
				for (int i = 0; i < var_PartC.size(); ++i) {
					iisphparticle& part = var_PartC[i];
				if(part.isboundary == true ){
				//if (((part.isboundary == true && part.ismovingboundary == true) || part.drawme == true) && part.pos.y < upperviualbord && part.pos.y > lowervisualbord ) {
					std::string xPos = std::to_string(part.pos.x);
					std::string yPos = std::to_string(part.pos.y);
					std::string zPos = std::to_string(part.pos.z);
					std::string isbound = std::to_string(part.ismovingboundary);
					if (part.ismovingboundary == false) {
						isbound = std::to_string(part.isfloatingboundary);
					}
					
					if (part.isobstacle == true) {
						isbound = std::to_string(1);
					}
					
					file_out << xPos << ";" << yPos << ";" << zPos << ";" << isbound << "\n";
				}
			}
			file_out.flush();
			file_out.close();
			std::cout << "frame" << std::to_string(numexport) << std::endl;
			animationtime = 0;
		}
		//debug--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		if (maxvel > overallmaxvel) {
			overallmaxvel = maxvel;
		}
		currentdens = 0;
		maxvel = 0;;
		//_____________________________________________________________________________________________________________________________________________________________________________________________________________

		//IISPH___________________________________________________________________________________________________________________________________________________________________________________________________
		if (iisph&& paussimul == false || (iisph && animationstep == 0)) {

			//IISPHalgorythm(ParticlesContainer, hashmap);
			//var_iiSPHalgorythm(ParticlesContainer, hashmap);
			var_iisph(var_PartC, hashmap);
			animationstep++;
			currenttimes.push_back(totaltimeanimated +deltaT);
			if (exportanimation) {
				animationtime += deltaT;
			}
		}
		//endIISPH___________________________________________________________________________________________________________________________________________________________________________________________________
		//SESPH___________________________________________________________________________________________________________________________________________________________________________________________________
		if ((ssph && paussimul == false  && TwoDsimul == false)|| (ssph && animationstep == 0 && TwoDsimul == false)) {
			ssphAlgo(var_PartC, hashmap);
			//SSPHalgorythm(ParticlesContainer, hashmap);
			animationstep++;
			currenttimes.push_back(totaltimeanimated +deltaT);
			if (exportanimation) {
				animationtime += deltaT;
			}
		}
		//endSESPH___________________________________________________________________________________________________________________________________________________________________________________________________
		//2DSESPH___________________________________________________________________________________________________________________________________________________________________________________________________
		if (TwoDsimul && paussimul == false || (twoDiisph && animationstep == 0)) {
			if (ssph) {
				ssphAlgotwoD(var_PartC, hashmap);
			}
			else {
				twoDiisph(var_PartC, hashmap);
			}
			
			//SSPHalgorythm(ParticlesContainer, hashmap);
			animationstep++;
			currenttimes.push_back(totaltimeanimated + deltaT);
			if (exportanimation) {
				animationtime += deltaT;
			}
		}
		//end2DSESPH___________________________________________________________________________________________________________________________________________________________________________________________________

		//debug--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		neighbouraverage = 0;
//#pragma omp parallel for
		for (int i = 0; i < var_MaxParticles; i++) {
			//var_PartC[i].drawme = false;
			if (var_PartC[i].isboundary) {
				var_PartC[i].r = 44;
				var_PartC[i].g = 2;
				var_PartC[i].r = 25;
				var_PartC[i].a = boundarya;
			}
			if (var_PartC[i].ismovingboundary || var_PartC[i].isfloatingboundary) {
				var_PartC[i].a = obstaclea;
			}
		}
		for (int i = 0; i < var_MaxParticles; i++) {
			overallminxpos = std::min(overallminxpos, var_PartC[i].pos.x);
			overallminypos = std::min(overallminypos, var_PartC[i].pos.y);
			overallminzpos = std::min(overallminzpos, var_PartC[i].pos.z);
			if (var_PartC[i].index == observing && observing < (var_MaxParticles + 1)) {
				observingposx = var_PartC[i].pos.x;
				observingposy = var_PartC[i].pos.y;
				observingposz = var_PartC[i].pos.z;

				var_PartC[i].r = 255;
				var_PartC[i].g = 0;
				var_PartC[i].b = 0;
				var_PartC[i].drawme = true;
				numNeighofobserving = 0;
				sumkernelderobserving = 0;
				sumkernelderobservingy = 0;
				sumkernelderobservingz = 0;
				for (std::tuple<int, glm::vec3, glm::vec3>& neigh : var_PartC[i].IdNSubKernelder) {
					if (var_PartC[std::get<0>(neigh)].index != observing) {
						var_PartC[std::get<0>(neigh)].r = 0;
						var_PartC[std::get<0>(neigh)].g = 255;
						var_PartC[std::get<0>(neigh)].b = 0;
						var_PartC[std::get<0>(neigh)].a = 255;
						var_PartC[std::get<0>(neigh)].drawme = true;
						numNeighofobserving += 1;
					}
					
					sumkernelderobserving += std::get<2>(neigh).x;
					sumkernelderobservingy += std::get<2>(neigh).y;
					sumkernelderobservingz += std::get<2>(neigh).z;
				}
				sumkernelobserving = 0;
				for (auto& neigh : var_PartC[i].Kernel) {
					sumkernelobserving += neigh;
				}
				densobserving = var_PartC[i].density;
				massobserving = var_PartC[i].m;
				presaobserving = glm::length(var_PartC[i].presA);
				nonpresaobserving = glm::length(var_PartC[i].nonpresA);
			}
			if (var_PartC[i].isboundary == false && var_PartC[i].isfloatingboundary == false) {
				for (std::tuple<int, glm::vec3, glm::vec3>& neigh : var_PartC[i].IdNSubKernelder) {
					neighbouraverage += 1;
				}
				if (glm::length(var_PartC[i].vel) > maxvel) {
					maxvel = glm::length(var_PartC[i].vel);
				}
				if (var_PartC[i].pos.y < minypos) {
					minypos = var_PartC[i].pos.y;
				}
				if (var_PartC[i].pos.y > maxypos) {
					maxypos = var_PartC[i].pos.y;
				}
			}
		}
		neighbouraverage = (neighbouraverage-var_fluidpart) / var_fluidpart;
		// Simulate all particles
		auto start_d = std::chrono::high_resolution_clock::now();
		int ParticlesCount = 0;

		for (iisphparticle& p : var_PartC) {
			if ((showallpart || (showboundandouter && (p.drawme == true || p.isboundary)) || (showouter&& (p.drawme == true) ))&& p.pos.y < upperviualbord && p.pos.y > lowervisualbord){
				g_particule_position_size_data[4 * ParticlesCount + 0] = p.pos.x * sizefac;
				g_particule_position_size_data[4 * ParticlesCount + 1] = p.pos.y * sizefac;
				g_particule_position_size_data[4 * ParticlesCount + 2] = p.pos.z * sizefac;
				g_particule_position_size_data[4 * ParticlesCount + 3] = p.size * sizefac*1;

				g_particule_color_data[4 * ParticlesCount + 0] = p.r;
				g_particule_color_data[4 * ParticlesCount + 1] = p.g;
				g_particule_color_data[4 * ParticlesCount + 2] = p.b;
				g_particule_color_data[4 * ParticlesCount + 3] = p.a;

				ParticlesCount++;
					
			}
		}
		animatedparticles = ParticlesCount;
		//_____________________________________________________________________________________________________________________________________________________________________________________________________________

				//printf("%d ",ParticlesCount);

				// Update the buffers that OpenGL uses for rendering.
				// There are much more sophisticated means to stream data from the CPU to the GPU, 
				// but this is outside the scope of this tutorial.
				// http://www.opengl.org/wiki/Buffer_Object_Streaming
		
		glBindBuffer(GL_ARRAY_BUFFER, particles_position_buffer);
		glBufferData(GL_ARRAY_BUFFER, var_MaxParticles * 4 * sizeof(GLfloat), NULL, GL_STREAM_DRAW); // Buffer orphaning, a common way to improve streaming perf. See above link for details.
		glBufferSubData(GL_ARRAY_BUFFER, 0, ParticlesCount * sizeof(GLfloat) * 4, g_particule_position_size_data);

		glBindBuffer(GL_ARRAY_BUFFER, particles_color_buffer);
		glBufferData(GL_ARRAY_BUFFER, var_MaxParticles * 4 * sizeof(GLubyte), NULL, GL_STREAM_DRAW); // Buffer orphaning, a common way to improve streaming perf. See above link for details.
		glBufferSubData(GL_ARRAY_BUFFER, 0, ParticlesCount * sizeof(GLubyte) * 4, g_particule_color_data);
		
		/*
		glBindBuffer(GL_ARRAY_BUFFER, particles_position_buffer);
		glBufferData(GL_ARRAY_BUFFER, ParticlesCount * 4 * sizeof(GLfloat), g_particule_position_size_data.data(), GL_STREAM_DRAW);

		glBindBuffer(GL_ARRAY_BUFFER, particles_color_buffer);
		glBufferData(GL_ARRAY_BUFFER, ParticlesCount * 4 * sizeof(GLubyte), g_particule_color_data.data(), GL_STREAM_DRAW);
		*/
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		// Use our shader
		glUseProgram(programID);

		// Bind our texture in Texture Unit 0
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, Texture);
		// Set our "myTextureSampler" sampler to use Texture Unit 0
		glUniform1i(TextureID, 0);

		// Same as the billboards tutorial
		glUniform3f(CameraRight_worldspace_ID, ViewMatrix[0][0], ViewMatrix[1][0], ViewMatrix[2][0]);
		glUniform3f(CameraUp_worldspace_ID, ViewMatrix[0][1], ViewMatrix[1][1], ViewMatrix[2][1]);

		glUniformMatrix4fv(ViewProjMatrixID, 1, GL_FALSE, &ViewProjectionMatrix[0][0]);

		// 1rst attribute buffer : vertices
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, billboard_vertex_buffer);
		glVertexAttribPointer(
			0,                  // attribute. No particular reason for 0, but must match the layout in the shader.
			3,                  // size
			GL_FLOAT,           // type
			GL_FALSE,           // normalized?
			0,                  // stride
			(void*)0            // array buffer offset
		);

		// 2nd attribute buffer : positions of particles' centers
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, particles_position_buffer);
		glVertexAttribPointer(
			1,                                // attribute. No particular reason for 1, but must match the layout in the shader.
			4,                                // size : x + y + z + size => 4
			GL_FLOAT,                         // type
			GL_FALSE,                         // normalized?
			0,                                // stride
			(void*)0                          // array buffer offset
		);

		// 3rd attribute buffer : particles' colors
		glEnableVertexAttribArray(2);
		glBindBuffer(GL_ARRAY_BUFFER, particles_color_buffer);
		glVertexAttribPointer(
			2,                                // attribute. No particular reason for 1, but must match the layout in the shader.
			4,                                // size : r + g + b + a => 4
			GL_UNSIGNED_BYTE,                 // type
			GL_TRUE,                          // normalized?    *** YES, this means that the unsigned char[4] will be accessible with a vec4 (floats) in the shader ***
			0,                                // stride
			(void*)0                          // array buffer offset
		);

		// These functions are specific to glDrawArrays*Instanced*.
		// The first parameter is the attribute buffer we're talking about.
		// The second parameter is the "rate at which generic vertex attributes advance when rendering multiple instances"
		// http://www.opengl.org/sdk/docs/man/xhtml/glVertexAttribDivisor.xml
		glVertexAttribDivisor(0, 0); // particles vertices : always reuse the same 4 vertices -> 0
		glVertexAttribDivisor(1, 1); // positions : one per quad (its center)                 -> 1
		glVertexAttribDivisor(2, 1); // color : one per quad                                  -> 1

		// Draw the particles !
		// This draws many times a small triangle_strip (which looks like a quad).
		// This is equivalent to :
		// for(i in ParticlesCount) : glDrawArrays(GL_TRIANGLE_STRIP, 0, 4), 
		// but faster.
		glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, ParticlesCount);

		glDisableVertexAttribArray(0);
		glDisableVertexAttribArray(1);
		glDisableVertexAttribArray(2);
		auto end_d = std::chrono::high_resolution_clock::now();
		auto duration_d = std::chrono::duration_cast<std::chrono::microseconds>(end_d - start_d);
		auto duration_alltime = std::chrono::duration_cast<std::chrono::microseconds>(end_d - start_alltime);
		if (oldanimstep != animationstep) {
			simstepcnter += 1;
			maxdenserrold = std::max(maxdenserrold, denserrold);
			overallmaxcfl = std::max(cfl, overallmaxcfl);
			drawingTime = duration_d.count();
			if (densitysnew.back()  > maxavgdensdeviation) {
				maxavgdensdeviation = densitysnew.back() ;
			}
			if (currentiter > maxiter) {
				maxiter = currentiter;
			}
			avgiter += currentiter;
			totalComputationTime = duration_alltime.count();
			
			totalComputationTime = totalComputationTime;
			if (totalComputationTime > maxcomptime) {
				maxcomptime = totalComputationTime;
			}
			avgcomptime += totalComputationTime;
			totaltimeanimated += deltaT;
			simcomptime = deltaT / (totalComputationTime / 1000000);
			avgsincomptime += simcomptime;
			fps = 1 / (totalComputationTime / 1000000);
			avgfps += fps;
			maxdenserralltime = std::max(maxdenserralltime, max_singledens);
			cfls.push_back(cfl);
			iterations.push_back(currentiter);
			densdiffes.push_back(glm::abs(densitysnew.back() - densitys.back()));
			maxdiffdens = std::max(densdiffes.back(), maxdiffdens);
			maxvels.push_back(maxvel);
		}
		if (resetvalues) {
			fps = 0;
			simstepcnter = 0;
			avgsincomptime = 0;
			avgfps = 0;
			maxdenserrold = 0;
			overallmaxcfl = 0;
			currenttimes.clear();
			densitysnew.clear();
			densitys.clear();
			cfls.clear();
			iterations.clear();
			densdiffes.clear();
			maxvels.clear();
			densitysnew.push_back(0);
			densitys.push_back(0);
			cfls.push_back(0);
			iterations.push_back(0);
			densdiffes.push_back(0);
			maxvels.push_back(0);
			maxavgdensdeviation = 0;
			currentiter = 0;
			maxiter = 0;
			avgiter = 0;
			totalComputationTime = 0;
			maxcomptime = 0;
			totaltimeanimated = 0;
			avgcomptime = 0;
			othercomputationTime = 0;
			simcomptime = 0;
			maxdenserralltime = 0;
			max_singledens = 0;
			minypos = 100000;
			maxypos = 0;
			overallmaxcfl = 0;
			overallmaxvel = 0;
			animationstep = 1;
			resetvalues = false;
			maxdiffdens = 0;
		}
		oldanimstep = animationstep;
		//imgui
		if (ImGui::Begin("Density Deviation Plot")) {
			// Plot für New Density Deviation
			ImGui::PlotLines("predicted Density Deviation", densitysnew.data(), densitysnew.size(), 0, nullptr, 0, 2.0f * denistyerrormax, ImVec2(0, 80));

			// Plot für Old Density Deviation
			ImGui::PushStyleColor(ImGuiCol_PlotLines, ImVec4(1, 0, 0, 1)); // Setze Farbe auf Rot
			ImGui::PlotLines("actual Density Deviation", densitys.data(), densitys.size(), 0, nullptr, 0, 2.0f * denistyerrormax, ImVec2(0, 80));
			ImGui::PopStyleColor();
			ImGui::PlotLines("difference Density Deviation",densdiffes.data(), densdiffes.size(), 0, nullptr, 0, maxdiffdens, ImVec2(0, 80));
			ImGui::PlotLines("cfl", cfls.data(), cfls.size(), 0, nullptr, 0, overallmaxcfl, ImVec2(0, 80));
			ImGui::PlotLines("solver iteration", iterations.data(), iterations.size(), 0, nullptr, 0, maxiter, ImVec2(0, 80));
			ImGui::PlotLines("maximum velocity", maxvels.data(), maxvels.size(), 0, nullptr, 0, overallmaxvel, ImVec2(0, 80));

		}
		ImGui::End();
		ImGui::Begin("schiebereglermachine");
		if (ImGui::Button("Step Simulation")) {
			StepSimulation(var_PartC, hashmap);
			paussimul = true;
		}
		ImGui::Checkbox("pause", &paussimul);
		ImGui::Checkbox("update Pos", &updateviewpos);
		ImGui::Text("fps: %.4f",fps);
		ImGui::Text("avg fps: %.4f",avgfps/simstepcnter);
		ImGui::Checkbox("reset stats", &resetvalues);
		ImGui::Text("estimated density deviation: %.3f %%", densitysnew.back() );
		ImGui::Text("density deviation: %.3f %%", densitys.back());
		ImGui::Text("max estimated density deviation: %.3f %%", maxavgdensdeviation);
		ImGui::Text("max density deviation: %.3f %%", maxdenserrold);
		ImGui::Text("used iterations: %i", currentiter);
		ImGui::Text("average comp time  %.1f ms", avgcomptime/animationstep);
		ImGui::Text("average iterations    %.1d", avgiter / animationstep);
		ImGui::Text("max comp time    %.1f ms",maxcomptime);
		ImGui::Text("max iterations    %.1d",maxiter);
		ImGui::Text("total simulation time    %.3f", totaltimeanimated);
		ImGui::Text("simtime per sec    %.3f", simcomptime);
		ImGui::Text("avg simtime per s    %.3f", avgsincomptime/simstepcnter);
		ImGui::Text("maximum density current: %.3f %%", max_singledens);
		ImGui::Text("maximum density ever: %.3f %%",maxdenserralltime);
		ImGui::Text("fluidparticles with a density lower then p0% .1d", var_fluidpart-numofp0high);
		ImGui::Text("fluidparticles unused for computation percent % .1d%%", 100 * (var_fluidpart - usemefordens.size()) / var_fluidpart);
		// Start a new ImGui table with 3 columns
		




		if (ImGui::CollapsingHeader("timerzz")) {
			if (ImGui::BeginTable("Time Metrics", 3, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg))
			{
				// Set column headers
				ImGui::TableSetupColumn("Timer");
				ImGui::TableSetupColumn("Time (ns)");
				ImGui::TableSetupColumn("Percentage (%)");
				ImGui::TableHeadersRow();

				// Row for Neighbour Search Time
				ImGui::TableNextRow();
				ImGui::TableSetColumnIndex(0);
				ImGui::Text("Neighbour Search Time");
				ImGui::TableSetColumnIndex(1);
				ImGui::Text("%.3f ns", neighSearchTime);
				ImGui::TableSetColumnIndex(2);
				ImGui::Text("%.1f%%", (neighSearchTime / (totalComputationTime )) * 100);

				// Row for Drawing Time
				ImGui::TableNextRow();
				ImGui::TableSetColumnIndex(0);
				ImGui::Text("Drawing Time");
				ImGui::TableSetColumnIndex(1);
				ImGui::Text("%.3f ns", drawingTime);
				ImGui::TableSetColumnIndex(2);
				ImGui::Text("%.1f%%", (drawingTime / (totalComputationTime )) * 100);

				// Row for Other Computation Time
				ImGui::TableNextRow();
				ImGui::TableSetColumnIndex(0);
				ImGui::Text("Kernel Computation Time");
				ImGui::TableSetColumnIndex(1);
				ImGui::Text("%.3f ns", kerneltime);
				ImGui::TableSetColumnIndex(2);
				ImGui::Text("%.1f%%", ((kerneltime) / (totalComputationTime )) * 100);
				// Row for Other Computation Time
				ImGui::TableNextRow();
				ImGui::TableSetColumnIndex(0);
				ImGui::Text("density Computation Time");
				ImGui::TableSetColumnIndex(1);
				ImGui::Text("%.3f ns", densitytime);
				ImGui::TableSetColumnIndex(2);
				ImGui::Text("%.1f%%", ((densitytime) / (totalComputationTime )) * 100);
				// Row for Other Computation Time
				ImGui::TableNextRow();
				ImGui::TableSetColumnIndex(0);
				ImGui::Text("non pressure acce Computation Time");
				ImGui::TableSetColumnIndex(1);
				ImGui::Text("%.3f ns", nonpresatime);
				ImGui::TableSetColumnIndex(2);
				ImGui::Text("%.1f%%", ((nonpresatime) / (totalComputationTime )) * 100);
				// Row for Other Computation Time
				ImGui::TableNextRow();
				ImGui::TableSetColumnIndex(0);
				ImGui::Text("pred vel Computation Time");
				ImGui::TableSetColumnIndex(1);
				ImGui::Text("%.3f ns", predveltime);
				ImGui::TableSetColumnIndex(2);
				ImGui::Text("%.1f%%", ((predveltime) / (totalComputationTime )) * 100);
				// Row for Other Computation Time
				ImGui::TableNextRow();
				ImGui::TableSetColumnIndex(0);
				ImGui::Text("sf Computation Time");
				ImGui::TableSetColumnIndex(1);
				ImGui::Text("%.3f ns", computeallsftime);
				ImGui::TableSetColumnIndex(2);
				ImGui::Text("%.1f%%", ((computeallsftime) / (totalComputationTime )) * 100);
				// Row for Other Computation Time
				ImGui::TableNextRow();
				ImGui::TableSetColumnIndex(0);
				ImGui::Text("aff Computation Time");
				ImGui::TableSetColumnIndex(1);
				ImGui::Text("%.3f ns", computeallafftime);
				ImGui::TableSetColumnIndex(2);
				ImGui::Text("%.1f%%", ((computeallafftime) / (totalComputationTime )) * 100);
				// Row for Other Computation Time
				ImGui::TableNextRow();
				ImGui::TableSetColumnIndex(0);
				ImGui::Text("loop Computation Time");
				ImGui::TableSetColumnIndex(1);
				ImGui::Text("%.3f ns", looptime);
				ImGui::TableSetColumnIndex(2);
				ImGui::Text("%.1f%%", ((looptime) / (totalComputationTime )) * 100);
				// Row for Total Computation Time
				ImGui::TableNextRow();
				ImGui::TableSetColumnIndex(0);
				ImGui::Text("ridgidbodytime");
				ImGui::TableSetColumnIndex(1);
				ImGui::Text("%.3f ns", othercomputationTime);
				ImGui::TableSetColumnIndex(2);
				ImGui::Text("%.1f%%", ((othercomputationTime) / (totalComputationTime)) * 100);

				// End the table
				ImGui::EndTable();
			}
		}
		if (ImGui::CollapsingHeader("integration setting")) {
			ImGui::Checkbox("iisph", &iisph);
			ImGui::Checkbox("ssph", &ssph);
			ImGui::Checkbox("uniform grid", &uniformgridneighsearch);
			ImGui::Checkbox("ssph pressure boundarys", &pressuredextrapolatebound);
			ImGui::Checkbox("2D iisph", &TwoDsimul);
		}
		if (ImGui::CollapsingHeader("scene")) {
			ImGui::Checkbox("add floating", &addfloating);
			ImGui::Checkbox("teslavalve", &realtesla);
			ImGui::Checkbox("obstacles open", &tesla);
			ImGui::Checkbox("obstacles closed", &teslaclosed);
			ImGui::Checkbox("small dambreak test", &smalldambreaktestscen);
			ImGui::Checkbox("big dambreak test", &dambreaktestscen);
			ImGui::Checkbox("moving boundary", &rotatinganimation);
			if (ImGui::CollapsingHeader("moving boundary settings")) {
				ImGui::Checkbox("single wall", &singlewall);
				ImGui::SliderInt("fluid in x direction", &var_nx, 1, var_oneway-2);
				ImGui::SliderInt("fluid in y direction", &var_ny, 1, var_oneway - 2);
				ImGui::SliderInt("fluid in z direction", &var_nz, 1, var_oneway - 2);
				ImGui::SliderInt("boundary cube", &var_oneway, 1, 55);
				ImGui::SliderInt("number of rotation parts", &num_rot_part, 50, 500);
				ImGui::SliderFloat("rotation speed", &baseRotationSpeed, 0, 0.8);
			}
			ImGui::Checkbox("2D dambreak", &TwoDDambreak);
			ImGui::Checkbox("2D water column", &TwoDwatercol);
			ImGui::Checkbox("small water column", &smallwatercol);
			if (ImGui::CollapsingHeader("small water column settings")) {
				ImGui::Checkbox("single wall", &singlewall);
				ImGui::InputInt("small watercolumn height", &watercolheight, 35, 210);
			}
			ImGui::Checkbox("water column", &watercol);
			if (ImGui::CollapsingHeader("water column settings")) {
				ImGui::Checkbox("single wall", &singlewall);
				ImGui::InputInt("watercolumn height", &watercolheight, 35, 210);
			}
			ImGui::Checkbox("Breaking dam", &start);
			if (ImGui::CollapsingHeader("breaking dam settings")) {
				ImGui::Checkbox("single wall", &singlewall);
				ImGui::SliderInt("fluid in x direction", &var_nx, 1, var_oneway - 2);
				ImGui::SliderInt("fluid in y direction", &var_ny, 1, var_oneway - 2);
				ImGui::SliderInt("fluid in z direction", &var_nz, 1, var_oneway - 2);
				ImGui::SliderInt("boundary cube", &var_oneway, 1, 55);
			}
			
		}
		if (ImGui::CollapsingHeader("simulation setting technical")) {
			ImGui::Checkbox("dont adapt dt to cfl", &fixeddt);
			ImGui::Checkbox("add pseudomass", &setboundmass);
			ImGui::Checkbox("add particlemass ", &setpartmass);
			ImGui::Checkbox("use absolut density deviation", &absinterrupt);
			ImGui::Checkbox("use the whole fluid for density deviation", &usewholefluidforcalc);
			ImGui::Checkbox("usewhole fluid for density division", &usewholefluidtodevide);
			ImGui::InputFloat("deltaTmax", &deltaTmax, 0.001f, 0.1f, "% .7f");
			ImGui::InputFloat("cflmax", &cfl_max, 0.1f, 1.f);
			ImGui::SliderFloat("deltaT", &deltaT, 0.00001f, 0.1f);
			ImGui::SliderFloat("k", &k, 10000, 10000000);
			ImGui::InputFloat("gamma1 density", &gammadens, 0.1, 2);
			ImGui::InputFloat("gamma2 pressure", &gammapres, 0.1, 2);
			ImGui::InputFloat("gamma3 bound", &gammabound, 0.1, 2);
			ImGui::InputFloat("jitterfactor ", &jitterfac, 0, 1);
			ImGui::InputFloat("rigidboundarymass", &massfac, 0.1, 2);
			ImGui::InputFloat("gamma particled", &gammapart, 0.1, 2);
			ImGui::InputFloat("omega ", &omega, 0.1, 1);
			ImGui::InputFloat("surfacetension ", &surfacetension, 0, 1);
			ImGui::SliderInt("maximum iterations ", &currentitermax, 1, 1000);
			ImGui::SliderInt("minimum iterations ", &currentitermin, 1, 1000);
			ImGui::SliderFloat("maximum density deviation ", &denistyerrormax,0.01,1);
			ImGui::SliderFloat("gravity ", &gravity, -19.81, 19.81);
			ImGui::InputFloat("boundary viscosity ", &visc_boundary, 0.001, 10, "% .8f");
			ImGui::InputFloat("particle viscosity ", &visc_fluid, 0.001, 10, "% .8f");

		}

		if (ImGui::CollapsingHeader("simulation setting visual")) {
			ImGui::SliderFloat("color multiplyer", &cloloroffset, 0, 50);
			ImGui::InputFloat("size", &sizefac, 0.1, 10000);
			ImGui::SliderInt("boundary alpha value", &boundarya, 0, 200);
			ImGui::SliderInt("obstacle alpha value", &obstaclea, 0, 200);
			ImGui::SliderInt("max neighbours ", &animationneighbourscount, 10, 60);
			ImGui::Checkbox("colorcode velocity", &colorvel);
			ImGui::Checkbox("colorcode density", &colordens);
			ImGui::Checkbox("colorcode pressure", &colorpres);
			ImGui::Checkbox("highlight extremefast", &highlightextremefast);
			ImGui::Checkbox("highlight underpopulated", &highlightunderpop);
			ImGui::Checkbox("Visualize all particles", &showallpart);
			ImGui::Checkbox("Visualize bound & outerparticles", &showboundandouter);
			ImGui::Checkbox("Visualize outerparticles", &showouter);
			ImGui::SliderFloat("visualizing upper boarder ", &upperviualbord, 50, 1000);
			ImGui::SliderFloat("visualizing lower boarder", &lowervisualbord, -100, -10);

		}
		if (ImGui::CollapsingHeader("overall stats")) {
			
			ImGui::Text("AllParticles: %i", var_MaxParticles);
			ImGui::Text("Fluidparticles: %i", var_fluidpart);
			ImGui::Text("Boundaryparticles: %i", var_MaxParticles- var_fluidpart);
			ImGui::Text("Animated particles: %i", animatedparticles);
			ImGui::Text("dt: %.8f", deltaT);
			ImGui::Text("cfl: %.8f", cfl);
			ImGui::Text("overall max cfl: %.8f", overallmaxcfl);
			
			ImGui::Text("average neighbours: %.3f", neighbouraverage);
			
			ImGui::Text("overall maximum velocity: %.1f", overallmaxvel);
			ImGui::Text("overall minimum ypos: %.1f", minypos/h);
			ImGui::Text("overall maximum ypos: %.1f", maxypos/h);
			ImGui::Text("current maximum velocity: %.1f", maxvel);
			ImGui::Text("current torque: %.1f", torque.length);
		}
		if (ImGui::CollapsingHeader("single stats")) {
			ImGui::Text("pos x: %.1f", observingposx/h);
			ImGui::Text("pos y % .1f", observingposy/h);
			ImGui::Text("pos z % .6f", observingposz/h);
			ImGui::Text("number of neighbours: %.1f", numNeighofobserving);
			ImGui::Text("density % .1f", densobserving);
			ImGui::Text("mass % .6f", massobserving);
			ImGui::Text("pressure accellatation % .6f", presaobserving);
			ImGui::Text("nonpressure accelaration % .6f", nonpresaobserving);
			ImGui::Text("sum kernel der x % .6f", sumkernelderobserving);
			ImGui::Text("sum kernel der y % .6f", sumkernelderobservingy);
			ImGui::Text("sum kernel der z % .6f", sumkernelderobservingz);
			ImGui::Text("sum kernel  % .6f", sumkernelobserving);
			ImGui::InputInt("osverving", &observing, 1, var_MaxParticles);
		}
		if (ImGui::CollapsingHeader("export")) {
			ImGui::Checkbox("export image", &exportimage);
			ImGui::Checkbox("export scene", &exportscene);
			ImGui::Checkbox("export densitys", &exportdens);
			ImGui::Checkbox("export animation", &exportanimation);
			ImGui::Checkbox("export convergence", &exportconvergence);
			ImGui::Checkbox("export cfl", &exportcfl);
			ImGui::InputText("ExportName", changename, 128);
		}
		if (ImGui::CollapsingHeader("import")) {
			ImGui::Checkbox("import scene", &importscene);
			ImGui::InputText("ImportName", changename, 128);
		}

		

		ImGui::End();
		ImGui::Render();
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
		//imgui end

		// Swap buffers
		glfwSwapBuffers(window);
		glfwPollEvents();
		glfwGetWindowSize(window, wWidth, wHeight);
		//GLFWmonitor* primaryMonitor = glfwGetPrimaryMonitor();
		//const GLFWvidmode* mode = glfwGetVideoMode(primaryMonitor);
		//glfwSetWindowMonitor(window, primaryMonitor, 0, 0, mode->width, mode->height, mode->refreshRate);

		glfwSetWindowSize(window, windowedWidth, windowedHeight);

	} // Check if the ESC key was pressed or the window was closed
	while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS &&
		glfwWindowShouldClose(window) == 0);


	delete[] g_particule_position_size_data;
	//g_particule_position_size_data.clear();
	//imgui
	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	//ImPlot::DestroyContext();
	ImGui::DestroyContext();
	//imgui end
	// Cleanup VBO and shader
	glDeleteBuffers(1, &particles_color_buffer);
	glDeleteBuffers(1, &particles_position_buffer);
	glDeleteBuffers(1, &billboard_vertex_buffer);
	glDeleteProgram(programID);
	glDeleteTextures(1, &Texture);
	glDeleteVertexArrays(1, &VertexArrayID);


	// Close OpenGL window and terminate GLFW
	glfwTerminate();

	return 0;
}