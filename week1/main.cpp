#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <chrono>
using namespace std::chrono;

#include "shader.h"
#include "camera.h"

#include "densityMap.h"
#include "data.h"
#include "time.h"

#define PI 3.141592653589

// If this is false, then the application will run in windowed mode
// If this is true, then the application will run in fullscren mode
#define FULLSCREEN 0

// Adjust these numbers depending on your monitor resolution
#if FULLSCREEN
#define SCR_WIDTH 1920
#define SCR_HEIGHT 1080
#else
#define SCR_WIDTH 800
#define SCR_HEIGHT 600
#endif

// Keyboard and mouse input functions
void cursorPosMovementCallback(GLFWwindow* window, double xpos, double ypos);
void cursorPosRotationCallback(GLFWwindow* window, double xpos, double ypos);
void processKeyboardInput(GLFWwindow* window);

// Updates the vertices on the graphics card
// -----
// This function is pretty slow right now (a few hundred milliseconds)
// because it writes several megabytes of data at once to the graphics card,
// but it will be optimized soon
void updateVertexBuffer(unsigned int& VBO, DensityMap& grid);

// Demo functions to show what the volume map looks like
void sphereDemo(DensityMap& grid);
void fanDemo(DensityMap& grid);
void realDemo(DensityMap& grid);

// Used in the mouse movement callback
double lastMouseX;
double lastMouseY;
bool firstMouse;

// Creating a Camera object
Camera cam;

int main() {
	// Window title
	std::string windowTitle = "Density Map";

    // Creating the density map
    int dim = 101;
    DensityMap grid(dim);

    // (Optional) Adds a fan-shaped arrangement of cells to the volume map
    realDemo(grid);

	// Initializing the OpenGL context
	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 4);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	// Creating the window object
	GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, windowTitle.c_str(), FULLSCREEN ? glfwGetPrimaryMonitor() : NULL, NULL);
	
	// If the window is not created (for any reason)
	if (window == NULL) {
		std::cout << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
		return -1;
	}

	// Setting callbacks
	glfwMakeContextCurrent(window);
	glfwSetCursorPosCallback(window, cursorPosMovementCallback);

	// Lock the cursor to the center of the window
	// and make it invisible
	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
	
	// Load the OpenGL functions from the graphics card
	if (!gladLoadGLLoader(GLADloadproc(glfwGetProcAddress))) {
		std::cout << "Failed to initialize GLAD" << std::endl;
		return -1;
	}

	// Creating the shaders for the cells in the cube
	// and for the lines of the border of the cube
    Shader cellShader("/home/yanwen/CML_CS/cml/week1/cells.vs", "/home/yanwen/CML_CS/cml/week1/cells.fs");
    Shader lineShader("/home/yanwen/CML_CS/cml/week1/lines.vs", "/home/yanwen/CML_CS/cml/week1/lines.fs");

	// Allows blending (translucent drawing)
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// Initializing mouse info
	lastMouseX = SCR_WIDTH / 2.0;
	lastMouseY = SCR_HEIGHT / 2.0;
	firstMouse = true;

	// Get the vertices from the volume map
	// in a form useful to OpenGL
	std::vector<float> cellPositions = grid.getVertices();
	std::vector<float> cellDensities = grid.getDensities();
	
	// Initializing the buffers storing the vertices
	// of the volume map on the graphics card
	unsigned int cellVAO, cellPositionVBO, cellDensityVBO;
	glGenBuffers(1, &cellPositionVBO);
	glGenBuffers(1, &cellDensityVBO);
	glGenVertexArrays(1, &cellVAO);

	glBindVertexArray(cellVAO);

	// Cell positions

	glBindBuffer(GL_ARRAY_BUFFER, cellPositionVBO);
	glBufferData(GL_ARRAY_BUFFER, cellPositions.size() * sizeof(float), cellPositions.data(), GL_STATIC_DRAW);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), 0);
	glEnableVertexAttribArray(0);

	// Data in each cell

	glBindBuffer(GL_ARRAY_BUFFER, cellDensityVBO);
	glBufferData(GL_ARRAY_BUFFER, cellDensities.size() * sizeof(float), cellDensities.data(), GL_STATIC_DRAW);
	
	glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, sizeof(float), 0);
	glEnableVertexAttribArray(1);

	// Array containing the coordinates of the vertices
	// of the white lines
	float lines[72] = {
		-5, -5, -5,
		 5, -5, -5,

		-5,  5, -5,
		 5,  5, -5,

		-5, -5,  5,
		 5, -5,  5,

		-5,  5,  5,
		 5,  5,  5,

		 // -----

		-5, -5, -5,
		-5,  5, -5,

		 5, -5, -5,
		 5,  5, -5,

		-5, -5,  5,
		-5,  5,  5,

		 5, -5,  5,
		 5,  5,  5,

		 // -----

		-5, -5, -5,
		-5, -5,  5,

		 5, -5, -5,
		 5, -5,  5,

		-5,  5, -5,
		-5,  5,  5,
		
		 5,  5, -5,
		 5,  5,  5,
	};

	// Initializing the buffer storing the vertices
	// of the white lines on the graphics card
	unsigned int lineVAO, lineVBO;
	glGenBuffers(1, &lineVBO);
	glGenVertexArrays(1, &lineVAO);

	glBindVertexArray(lineVAO);

	glBindBuffer(GL_ARRAY_BUFFER, lineVBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(lines), lines, GL_STATIC_DRAW);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), 0);
	glEnableVertexAttribArray(0);

	// Main event loop
	auto t1 = high_resolution_clock::now();
	int count = 0;
	while (!glfwWindowShouldClose(window)) {
	    ++count;
        if (count == 333){
            auto t2 = high_resolution_clock::now();
            auto duration = duration_cast<microseconds>(t2 - t1);
            std::cout << "time of 333 loops in ms is " << duration.count() << std::endl;
        }
		double currentFrame = glfwGetTime();
		cam.deltaTime = currentFrame - cam.lastFrame;
		cam.lastFrame = currentFrame;

		// Self-explanatory
		processKeyboardInput(window);

		// Clears the screen and fills it a dark grey color
		glClearColor(0.1, 0.1, 0.1, 1.0);
		glClear(GL_COLOR_BUFFER_BIT);

		// Creating matrices to transform the vertices into NDC (screen) coordinates
		// between -1 and 1 that OpenGL can use
		glm::dmat4 projection = glm::perspective(glm::radians(cam.fov), double(SCR_WIDTH) / SCR_HEIGHT, 0.01, 500.0);
		glm::dmat4 camView = cam.getViewMatrix();

		int dim = grid.getDim();
		glm::dmat4 model = glm::scale(glm::dmat4{}, glm::dvec3(10.0 / (dim - 1), 10.0 / (dim - 1), 10.0 / (dim - 1)));
		model = glm::translate(model, glm::dvec3(-(dim - 1) / 2.0, -(dim - 1) / 2.0, -(dim - 1) / 2.0));

		// Drawing the volume map
		cellShader.use();
		cellShader.setMat4("projection", projection);
		cellShader.setMat4("view", camView);
		cellShader.setMat4("model", model);

		glBindVertexArray(cellVAO);
		glDrawArrays(GL_TRIANGLES, 0, cellDensities.size());

		// Drawing the white lines
		lineShader.use();
		lineShader.setMat4("projection", projection);
		lineShader.setMat4("view", camView);
		lineShader.setMat4("model", glm::dmat4{});

		glBindVertexArray(lineVAO);
		glDrawArrays(GL_LINES, 0, 24);

		cam.prevPos = cam.position;

		// Update the screen
		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	// GLFW cleanup
	glfwTerminate();

	// Stops the window from closing immediately
	system("pause");
}

void processKeyboardInput(GLFWwindow *window) {
	// If shift is held down, then the camera moves faster
	bool sprinting = glfwGetKey(window, GLFW_KEY_LEFT_SHIFT);

	// WASD + Q and E movement
	if (glfwGetKey(window, GLFW_KEY_ESCAPE))
		glfwSetWindowShouldClose(window, true);
	if (glfwGetKey(window, GLFW_KEY_W))
		cam.processKeyboard(FORWARD, sprinting);
	if (glfwGetKey(window, GLFW_KEY_S))
		cam.processKeyboard(BACKWARD, sprinting);
	if (glfwGetKey(window, GLFW_KEY_A))
		cam.processKeyboard(LEFT, sprinting);
	if (glfwGetKey(window, GLFW_KEY_D))
		cam.processKeyboard(RIGHT, sprinting);
	if (glfwGetKey(window, GLFW_KEY_Q))
		cam.processKeyboard(DOWN, sprinting);
	if (glfwGetKey(window, GLFW_KEY_E))
		cam.processKeyboard(UP, sprinting);

	// Hold C to zoom in
	if (glfwGetKey(window, GLFW_KEY_C)) {
		cam.fov = 30;
	}
	else {
		cam.fov = 70;
	}
}

void cursorPosMovementCallback(GLFWwindow* window, double xpos, double ypos) {
	// Ensures that the viewer faces forward on startup
	if (firstMouse) {
		lastMouseX = xpos;
		lastMouseY = ypos;
		firstMouse = false;
	}

	// Updating the camera angle
	double xoffset = xpos - lastMouseX;
	double yoffset = lastMouseY - ypos;
	lastMouseX = xpos;
	lastMouseY = ypos;

	cam.processMouseMovement(xoffset, yoffset);
}

void sphereDemo(DensityMap& grid) {
	// Adds a sphere to the center of the volume map

	int dim = grid.getDim();

	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			for (int k = 0; k < dim; k++) {
				float xd = i - ((dim - 1) / 2.0);
				float yd = j - ((dim - 1) / 2.0);
				float zd = k - ((dim - 1) / 2.0);

				float mxd = (dim - 1) / 2.0;
				float myd = (dim - 1) / 2.0;
				float mzd = (dim - 1) / 2.0;

				float distance = sqrt(xd * xd + yd * yd + zd * zd);
				float maxDistance = sqrt(mxd * mxd + myd * myd + mzd * mzd);
				float shade = (maxDistance - distance) / maxDistance;
				shade = shade * shade;

				grid.cells[i][j][k] = shade;
			}
		}
	}
}

void fanDemo(DensityMap& grid) {
	// Adds a fan shape to the volume map
	// using the DensityMap::addLine() function

	glm::vec3 vertex = { 0.5, 0.5, 0.5 };

	float a1 = 1;
	float a2 = 1;

	float r = 0.3;

	for (; a2 <= 3; a2 += 0.01) {
		float x = r * sin(a1) * cos(a2);
		float y = r * sin(a1) * sin(a2);
		float z = r * cos(a1);

		std::vector<float> vals;

		for (int i = 0; i < 1000; i++) {
			vals.push_back(1);
		}

		//grid.addLine(vertex, vertex + glm::vec3(x, y, z), vals);
		grid.addLineSmoothed(vertex, vertex + glm::vec3(x, y, z), vals, 3);
	}
}

void realDemo(DensityMap& grid){
    //read the data from current red pitaya 2d data.
    std::vector<int> marker_locations;
    std::vector<scan_data_struct> scan_data;
    std::vector<screen_data_struct> screen_data;

    std::ifstream inFile("/home/yanwen/Downloads/clear_1.txt", std::ios::in | std::ios::binary);
    if (!inFile){
        printf("Failed to open file.\n");
        //return -1;
    }
    /* convert file to bytes vector */
    /* DO NOT USE ISTREAM_ITERATOR*/
    std::vector<unsigned char> file_bytes(
            (std::istreambuf_iterator<char>(inFile)),
            (std::istreambuf_iterator<char>()));
    for(i = 0; i < 20; i++){
        printf("%02X ", file_bytes.at(i));
    }
    printf("\n");
    /* find all marker locations */
    marker_locations = find_marker(file_bytes);
    /* convert file bytes to data struct */
    file_to_data(file_bytes, marker_locations, scan_data);
    /* convert data to vertex on screen */
    data_to_pixel(scan_data, screen_data);
    printf("find the screen_data\n");

    //fill the cell
    double maxx = 0, maxy = 0, maxz = 0, minx = 0, miny = 0, minz = 0;

    for (auto s: screen_data) {
        maxx = std::max(s.X, maxx);
        maxy = std::max(s.Y, maxy);
        maxz = std::max(s.Z, maxz);
        minx = std::min(s.X, minx);
        miny = std::min(s.Y, miny);
        minz = std::min(s.Z, minz);
    }

    /*
    for (i = 0; i < screen_data.size(); ++i)
    {
        maxx = screen_data.at(i).X > maxx? screen_data.at(i).X : maxx;
        maxy = screen_data.at(i).Y > maxy? screen_data.at(i).Y : maxy;
        maxz = screen_data.at(i).Z > maxz? screen_data.at(i).Z : maxz;
        minx = screen_data.at(i).X < minx? screen_data.at(i).X : minx;
        miny = screen_data.at(i).Y < miny? screen_data.at(i).Y : miny;
        minz = screen_data.at(i).Z < minz? screen_data.at(i).Z : minz;
    }
    */
    printf("find the maxes and mins\n");
    int ddim = grid.getDim();
    for (auto s: screen_data)
    {
        int tx = (int) ((s.X-minx) / (maxx-minx) * ddim);
        int ty = (int) ( (s.Y - miny) / (maxx-minx) * ddim);
        //int tz = (int) ((s.Z - minz) / (maxx-minx) * ddim);
        int tz = ddim/2;
        if (tx < 0 || tx >= ddim || ty < 0 || ty >= ddim || tz < 0 || tz >= ddim)
            continue;
        grid.cells[tx][ty][tz] = (s.I + 0.5);
    }
}

void updateVertexBuffer(unsigned int& VBO, DensityMap& grid) {
	// Gets the vertices from the density map
	std::vector<float> densities = grid.getDensities();

	// Writes the vertices to the vertex buffer on the graphics card
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferSubData(GL_ARRAY_BUFFER, 0, densities.size() * sizeof(float), densities.data());
}

double map_range_to_range(double input_start, double input_end, double output_start, double output_end, double input){
    return output_start + ((output_end - output_start) / (input_end - input_start)) * (input - input_start);
}

double convert_angle_2d_probe(double angle){
    double left = 353; double right = 176; double top = 85; double bottom = 262;
    if ((angle >= right) && (angle < bottom)){
        return map_range_to_range(right, bottom, -120, -90, angle);
    }
    else if ((angle >= top) && (angle < right)){
        return map_range_to_range(top, right, -90, -120, angle);
    }
    else if ((angle >= left) || (angle < top)){
        if ((angle <= top)){angle = angle+360;}
        return map_range_to_range(left, top+360, -60, -90, angle);
    }
    else if ((angle >= bottom) && (angle <= left)){
        return map_range_to_range(left, bottom, -60, -90, angle);
    }
}


void data_to_pixel(std::vector<scan_data_struct> _scan_data, std::vector<screen_data_struct> & _screen_data){
    //printf("%d\n", (int)_scan_data.size());
    for (i = 0; i < (int)_scan_data.size(); i++){
        double angle = _scan_data.at(i).encoder * 360.0 / 4096.0;
        angle = convert_angle_2d_probe(angle);
        /* find min and max */
        for (j = 0; j < buffer_length; j++){
            adc_max = std::max(adc_max, _scan_data.at(i).buffer[j]);
            adc_min = std::min(adc_min, _scan_data.at(i).buffer[j]);
        }
        /* normalize on the go */
        for (j = 0; j < buffer_length; j++){
            intensity = ((double)_scan_data.at(i).buffer[j] - adc_min)/(adc_max-adc_min);
            screen_data_struct temp_data = {(j+1) * Cos(angle), (j+1) * Sin(angle), 0, intensity};
            _screen_data.push_back(temp_data);
        }
        adc_max = 0; adc_min = 0;
    }
}

void file_to_data(std::vector<unsigned char> _file_bytes, std::vector<int> _marker_locations, std::vector<scan_data_struct> & _scan_data){
    for (i = 0; i < (int)_marker_locations.size()-1; i++){
        marker_index = _marker_locations.at(i);
        marker_index_next = _marker_locations.at(i+1);
        /* time stamp */
        for (j = 0; j < (int)sizeof(time_stamp_char); j++){
            time_stamp_char[j] = _file_bytes.at(marker_index + sizeof(marker) + j);
        }
        std::memcpy(&time_stamp, time_stamp_char, sizeof(time_stamp));
        time_stamp = changed_endian_4Bytes(time_stamp);
        /* probe type char */
        probe_type_char = _file_bytes.at(marker_index + sizeof(marker) + sizeof(time_stamp_char));
        /* encoder */
        for (j = 0; j < (int) sizeof(encoder_char); j++){
            encoder_char[j] = _file_bytes.at(marker_index + sizeof(marker) + sizeof(time_stamp_char) + sizeof(probe_type_char) + j);
        }
        std::memcpy(&encoder, encoder_char, sizeof(encoder));
        encoder = changed_endian_2Bytes(encoder);
        /* adc */
        /* determine the length of buffer */
        buffer_length = (int)(_marker_locations.at(i+1) - _marker_locations.at(i) - sizeof(marker) - sizeof(time_stamp_char) -
                              sizeof(probe_type_char) - sizeof(encoder_char) - sizeof(crc_char))/2;
        for (j = 0; j < buffer_length; j++){
            for (k = 0; k < (int)sizeof(adc_temp); k++){
                adc_temp[k] = _file_bytes.at(marker_index + sizeof(marker) + sizeof(time_stamp_char) + sizeof(probe_type_char) + sizeof(encoder_char) + j * 2 + k);
                adc_char[2*j+k] = adc_temp[k];
            }
            std::memcpy(&adc, adc_temp, sizeof(adc));
            adc = changed_endian_2Bytes(adc);
            buffer[j] = adc;
        }
        /* checksum */
        for (j = 0; j < (int)sizeof(crc_char); j++){
            crc_char[j] = _file_bytes.at(marker_index_next-(int)sizeof(crc_char)+j);
        }
        /* calculate crc locally */
        memcpy(crc_input, time_stamp_char, sizeof(time_stamp_char));
        memcpy(crc_input+sizeof(time_stamp_char), &probe_type_char, sizeof(probe_type_char));
        memcpy(crc_input+sizeof(time_stamp_char)+sizeof(probe_type_char), encoder_char, sizeof(encoder_char));
        memcpy(crc_input+sizeof(time_stamp_char)+sizeof(probe_type_char)+sizeof(encoder_char), adc_char, sizeof(adc_char));
        crc_result = crc32c(0, crc_input, sizeof(crc_input));
        crc_result = changed_endian_4Bytes(crc_result);
        memcpy(crc_result_char, (unsigned char *)&crc_result, sizeof (crc_result));
        /* if two crc matches */
        if (compare_crc(crc_char, crc_result_char, sizeof(crc_char))){
            scan_data_struct temp_struct;
            temp_struct.time_stamp = time_stamp;
            temp_struct.encoder = encoder;
            /* normalize on the go */
            for (j = 0; j < buffer_length; j++) {
                temp_struct.buffer[j] = buffer[j];
                //printf("Intensity:%f\n", temp_struct.buffer[j]);
            }
            _scan_data.push_back(temp_struct);
        }
    }
}

int compare_crc(unsigned char a[], unsigned char b[], size_t len){
    int ii;
    for (ii = 0; ii < (int)len; ii++){
        if (a[ii] != b[ii]){
            return 0;
        }
    }
    return 1;
}


std::vector<int> find_marker(std::vector<unsigned char> _file_bytes){
    std::vector<int> _marker_locations;
    for (i = 0; i < (int)_file_bytes.size(); i++){
        if (_file_bytes.at(i) == marker[0]){
            for (j = 0; j < (int)sizeof(marker); j++) {
                if ((i + j) < _file_bytes.size()){
                    if (_file_bytes.at(i + j) != marker[j]) {
                        marker_flag = false;
                        break;
                    } else {
                        marker_flag = true;
                    }
                }
            }
            if (marker_flag) {
                _marker_locations.push_back(i);
                i += (int)sizeof(marker);
            }
        }
        marker_flag = false;
    }
    return _marker_locations;
}

unsigned long changed_endian_4Bytes(unsigned long num){
    int byte0, byte1, byte2, byte3;
    byte0 = (num & 0x000000FF) >> 0 ;
    byte1 = (num & 0x0000FF00) >> 8 ;
    byte2 = (num & 0x00FF0000) >> 16 ;
    byte3 = (num & 0xFF000000) >> 24 ;
    return((byte0 << 24) | (byte1 << 16) | (byte2 << 8) | (byte3 << 0));
}

int16_t changed_endian_2Bytes(int16_t value){
    return ((value >> 8) & 0x00ff) | ((value & 0x00ff) << 8);
}

uint32_t crc32c(uint32_t crc, const unsigned char *buf, size_t len)
{
    int k;

    crc = ~crc;
    while (len--) {
        crc ^= *buf++;
        for (k = 0; k < 8; k++)
            crc = crc & 1 ? (crc >> 1) ^ 0xedb88320 : crc >> 1;
    }
    return ~crc;
}
