#pragma once
#include <random>
#include "ofMain.h"

class ofApp : public ofBaseApp{
	public:
		void setup();
		void update();
		void draw();
		
		void keyPressed(int key);
		void keyReleased(int key);
		void mouseMoved(int x, int y);
		void mouseDragged(int x, int y, int button);
		void mousePressed(int x, int y, int button);
		void mouseReleased(int x, int y, int button);
		void windowResized(int w, int h);
		void dragEvent(ofDragInfo dragInfo);
		void gotMessage(ofMessage msg);
    
        ofShader shader;
        ofFbo fbo;
        ofFbo fbo2;
        bool pass = false;
        bool diffuse = false;
        bool stop = false;
        bool mouseClick = false;
        float h = 0.01;
        glm::vec3 camPos;
        glm::vec3 mcamPos;
        glm::vec3 scamPos;
        glm::vec2 mosPos;
        glm::vec3 hiscamPos[256];
        uint32_t iters=0;
        std::default_random_engine generator;
        std::normal_distribution<double> distribution{0.0,1.0};
        ofVideoGrabber vidGrabber;
        ofFbo camfbo;
        ofPixels videoInverted;
        ofTexture videoTexture;
        int camWidth;
        int camHeight;
};
