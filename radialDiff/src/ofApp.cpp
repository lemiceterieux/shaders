#include "ofApp.h"

//--------------------------------------------------------------
void ofApp::setup(){
    
#ifdef TARGET_OPENGLES
	shader.load("shadersES2/shader");
#else
	if(ofIsGLProgrammableRenderer()){
		shader.load("shadersGL3/shader");
	}else{
		shader.load("shadersGL2/shader");
	}
#endif

    ofDisableArbTex();
	fbo.allocate(1280,800, GL_RGBA);
    fbo.begin();
    ofClear(255,255,255, 0);
    fbo.end();
    ofDisableArbTex();
	fbo2.allocate(1280,800,GL_RGBA);
    fbo2.begin();
    ofClear(255,255,255, 0);
    fbo2.end();
}

//--------------------------------------------------------------
void ofApp::update(){
    //
}

//--------------------------------------------------------------
void ofApp::draw(){
    ofSetColor(255);
    

    fbo.begin();
    shader.begin();
    shader.setUniform1f("time", ofGetElapsedTimef());
    shader.setUniform1f("h",h);
    shader.setUniform2f("resolution",1280, 800);
    shader.setUniform1i("diffuse",diffuse);
    shader.setUniformTexture("prevFrame", fbo2.getTexture(), 1);
    ofDrawRectangle(0, 0, 1280, 800);

    shader.end();
    fbo.end();
    if(pass){
        fbo.draw(0,0);
    }

    fbo2.begin();
    shader.begin();
    shader.setUniform1f("time", ofGetElapsedTimef());
    shader.setUniform1f("h",h);
    shader.setUniform2f("resolution",1280, 800);
    shader.setUniform1i("diffuse",diffuse);
    shader.setUniformTexture("prevFrame", fbo.getTexture(), 1);
    ofDrawRectangle(0, 0, 1280, 800);
   
    shader.end();
    fbo2.end();
    if(!pass){
        fbo2.draw(0,0);
    }
    pass = !pass;

}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){
    if(key == ' ') {
        diffuse = !diffuse;
    }
    if(key == OF_KEY_UP) {
        h *=2;
    }
    if(key == OF_KEY_DOWN) {
        h /=2;
    }

}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y){

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}
