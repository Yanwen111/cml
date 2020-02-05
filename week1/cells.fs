// FRAGMENT SHADER

#version 440 core

out vec4 FragColor;

in float fShade;

void main() {
	float shade = pow(fShade, 4.0) * abs(fShade);
	shade = clamp(shade, 0.0025, 1.0);
	FragColor = vec4(1.0, 1.0, 1.0, shade);
}