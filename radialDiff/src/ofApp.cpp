#include "ofApp.h"
#define CAPTURE_HEIGHT  480
#define CAPTURE_WIDTH   640

glm::vec3 integrate(glm::vec3 data){
    float sigma = 10;
    float rho = 28;
    float beta = 8/3;
    float diff[3];
    glm::vec3 ret(0,0,0);
    diff[0] = sigma*(data[1] - data[0]);
    diff[1] = data[0]*(rho-data[2])-data[1];
    diff[2] = data[0]*(data[1])-beta*data[2];
    ret.x = data.x +0.01*diff[0];
    ret.y = data.y +0.01*diff[1];
    ret.z = data.z +0.01*diff[2];
    return ret;
}

//--------------------------------------------------------------
void ofApp::setup(){
    camWidth = 640;  // try to grab at this size.
    camHeight = 320;
    vidGrabber.setDeviceID(0);
    vidGrabber.setDesiredFrameRate(30);
    vidGrabber.initGrabber(camWidth, camHeight);    
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
    ofDisableArbTex();
	camfbo.allocate(CAPTURE_WIDTH,CAPTURE_HEIGHT,GL_RGBA);
    camfbo.begin();
    ofClear(255,255,255, 0);
    camfbo.end();


    camPos = glm::vec3(1,1,1);
    mcamPos = glm::vec3(0,0,0);
    scamPos = glm::vec3(0,0,0);
//    hiscamPos = glm::vec3[256];
}

//--------------------------------------------------------------
void ofApp::update(){
    //
    glm::vec3 tempcamPos = integrate(camPos);
    hiscamPos[iters++%256] = tempcamPos;
    mcamPos = glm::vec3(0,0,0);
    for (int i = 0; i < 256; i++){
        mcamPos += hiscamPos[i];
        scamPos += hiscamPos[i]*hiscamPos[i];
    }
    mcamPos /= 256;
    scamPos /= 256;
    camPos = tempcamPos;
}

//--------------------------------------------------------------
void ofApp::draw(){
    ofSetColor(255);
    
    if(mouseClick){
        mosPos.x = mouseX;
        mosPos.y = mouseY;
        std::cout<<mosPos.x << std::endl << mosPos.y << std::endl;
    }

    vidGrabber.update();
    camfbo.begin();
    vidGrabber.draw(0, 0, CAPTURE_WIDTH, CAPTURE_HEIGHT);  
    camfbo.end();
    glm::vec3 heat(distribution(generator),distribution(generator),distribution(generator));
    fbo.begin();
    shader.begin();
    shader.setUniform1f("time", ofGetElapsedTimef());
    shader.setUniform1f("h",h);
    shader.setUniform2f("resolution",1280,800);
    shader.setUniform3f("heat", heat);
    shader.setUniform3f("camPos",camPos);
    shader.setUniform3f("mcamPos",mcamPos);
    shader.setUniform3f("scamPos",scamPos);
    shader.setUniform2f("mouse",mosPos);
    shader.setUniform1i("diffuse",diffuse);
    shader.setUniform1i("stop",stop);
    shader.setUniformTexture("prevFrame", fbo2.getTexture(), 0);
    shader.setUniformTexture("video", camfbo.getTexture(),1);
    ofDrawRectangle(0, 0, 1280,800);

    shader.end();
    fbo.end();
    if(pass){
        fbo.draw(0,0);
    }

    fbo2.begin();
    shader.begin();
    shader.setUniform1f("time", ofGetElapsedTimef());
    shader.setUniform1f("h",h);
    shader.setUniform2f("resolution",1280,800);
    shader.setUniform3f("heat", heat);
    shader.setUniform3f("camPos",camPos);
    shader.setUniform3f("mcamPos",mcamPos);
    shader.setUniform3f("scamPos",scamPos);
    shader.setUniform2f("mouse",mosPos);
    shader.setUniform1i("diffuse",diffuse);
    shader.setUniform1i("stop",stop);
    shader.setUniformTexture("prevFrame", fbo.getTexture(), 0);
    shader.setUniformTexture("video", camfbo.getTexture(),1);
    ofDrawRectangle(0, 0, 1280,800);
   
    shader.end();
    fbo2.end();
    if(!pass){
        fbo2.draw(0,0);
    }
    pass = !pass;
//    ofDrawBitmapString("Drag spacebar for diffusion", 15,15);
//    ofDrawBitmapString("Press UP or DOWN to change integration time", 15, 30);
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){
    if(key == ' ') {
        diffuse = !diffuse;
    }
    if(key == 'n') {
        stop = !stop;
    }
    if(key == OF_KEY_UP) {
        std::cout << "HELLO" << std::endl;
        h *=2;
    }
    if(key == OF_KEY_DOWN) {
        h /=2;
    }
    if(key == 'c'){
        shader.load("shadersGL3/shader");
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
    mouseClick= true;
}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){
    mouseClick= false;
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
