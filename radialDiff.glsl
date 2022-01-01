#version 150
#ifdef GL_ES
    precision mediump float;
#endif
uniform float time;
uniform vec2 resolution;
uniform vec2 mouse;

uniform sampler2D texture0;
uniform sampler2D texture1;
uniform sampler2D texture2;
uniform sampler2D texture3;
uniform sampler2D prevFrame;
uniform sampler2D prevPass;

#define PI 3.141592653

// Vertex data to manipulate
in VertexData
{
    vec4 v_position;
    vec3 v_normal;
    vec2 v_texcoord;
} inData;

// Struct for holding SDF return and color
struct sdCol{
    float sd;
    vec3 col;
};

// Matrix multiplation can come in handy someday
vec3 matMult(mat3 Map, vec3 vec){
    vec3 temp = vec;
    temp.x = dot(Map[0],temp);
    temp.y = dot(Map[1],temp);
    temp.z = dot(Map[2],temp);
    return temp;
}

// pseudorandom...
float rand(vec2 co){
    return fract(sin(dot(co, vec2(12.9898, 78.233))) * 43758.5453);
}

out vec4 fragColor;

// Move sphere 1 along cylindrical coordinates; Sphere 2 over cartesian
float RP = exp(cos(time)-1);
float thetaP = mod(time,2*3.14);
float phiP = mod(time-2*3.14,2*3.14);
vec3 pos = vec3(RP*sin(thetaP)*cos(phiP),RP*sin(thetaP)*sin(phiP),1.3-.1*.1*sin(time));
vec3 pos2 = vec3(1*cos(time),1*sin(cos(sin(time)))+.1,.1*cos(time)+1);

// Ground object
sdCol sdGround(vec3 p){
    sdCol sdG;
    sdG.sd = 0;
    // Mirror over x and y dimensions
    // Threshold on Z dimension
    vec3 sdd = vec3(cos((-(p.x))),(cos(1*(p.y))),((exp(-p.z))));
    sdG.sd = 2*(sdd.x + sdd.y + sdd.z);
    sdG.col = vec3(1*sin(p.x)*sin(time)+1.1,sin(time)*(p.y)+1.1,exp(-p.z));
    return sdG;
}

// Sphere object
sdCol sdSpheres(vec3 p, vec3 cent, float rad, vec3 col){
    sdCol sdS;
    sdS.sd = length(p -cent) - rad;
    sdS.col = col;
    return sdS;
}


// Rectangular prism object maybe
sdCol sdRectPrism(vec3 p, vec3 cent, vec3 dim, vec3 col){
    sdCol sdS;
    dim = dim;
    vec3 surfDist = (abs(p/cent) - dim);
    sdS.sd = length(max(surfDist,0) + min(max(surfDist.x, max(surfDist.y,surfDist.z)),0));
    
    sdS.col = col;
    return sdS;
}

// Minimum of an SDF and returning also the color of that object in a struct
sdCol minWithCol(sdCol obj1, sdCol obj2){
    if (obj1.sd < obj2.sd){
        return obj1;
    }else{
        return obj2;
    }
}

// Some AND operation on objects' SDF
sdCol addShad(sdCol obj1,sdCol obj2){
    sdCol temp = obj1;
    temp.sd = obj1.sd*obj2.sd;
    temp.col = obj1.col/(obj2.col+1*pow(10,-5));
    return temp;
}

// Our scene is two spheres and our ground plane
sdCol mapofWorld(vec3 p){

    sdCol ground = sdGround(p);
    sdCol Sphere1 = sdSpheres(p, pos, .3, 4*vec3(1*p.z,2*p.y,0.1*cos(time*p.x)));    
    sdCol Sphere2 = sdSpheres(p, pos2, .2, 3*vec3(cos(0)*cos(sin(p.z))*.1,.4*sin(p.z),cos(0)));
    sdCol Rect = sdRectPrism(p, pos2, vec3(.2), 3*vec3(cos(0)*cos(sin(p.z))*.1,.4*sin(p.z),cos(0)));    
    sdCol tempCol = minWithCol(Sphere1,Sphere2);
    tempCol = minWithCol(tempCol, ground);
    return tempCol;
}


const float MAXIMUM_TRACE_DISTANCE = 1000.0;

// Calculating normal vector to our object  by calculating the gradient of the signed distance
// function with respect to the ray position. 
vec3 calcNorm(vec3 p){
    const vec3 small_step = vec3(0.001, 0.00, 0.00);
    float gradient_x = mapofWorld(p + small_step.xyy).sd - mapofWorld(p - small_step.xyy).sd;
    float gradient_y = mapofWorld(p + small_step.yxy).sd - mapofWorld(p - small_step.yxy).sd;
    float gradient_z = mapofWorld(p + small_step.yyx).sd - mapofWorld(p - small_step.yyx).sd;

    vec3 normal = vec3(gradient_x, gradient_y, gradient_z);

    return normalize(normal);
}

// Adaptive distance ray marching
sdCol ray(vec3 r, vec3 rd){
    vec3 ro = r;  
    float dO = 0;
    float tol = 1e-5;
    float stepSize = .01;
    for (int i = 0; i < 16; i++){
        ro += ro + rd*(dO);
        sdCol dist = mapofWorld(ro);
        dO += dist.sd*.75;        
        if (dist.sd < tol){
            dist.sd = dO;
            return dist;
        }
        if (dO >= MAXIMUM_TRACE_DISTANCE){
            break;
        }
    }
    sdCol zero;
    zero.sd = 0;
    zero.col = vec3(0);
    return zero;
}

// Get shadow by ray marching from the normal of the hit object and see when we
// hit something. if we hit something, fill in shadow, otherwise scale  by
// distance to our object for soft shadow
float rayShad(vec3 r, vec3 rd){
    vec3 ro = r;  
    float dO = 0;
    float tol = 1e-4;
    float tmint = .1;
    float tmaxt = 1;
    float res = 1;
    int k = 32;
    for (float t = tmint; t < tmaxt;){
        ro += ro + rd*(t);
        sdCol dist = mapofWorld(ro);
        t += dist.sd;
        res = min(res,k*dist.sd/t);
        if (dist.sd < tol){
            dist.sd = dO;
            return 0;
        }
    }
    return res;//back/10;
}

float getLighting(vec3 p){
    vec3 norm = calcNorm(p);
    vec3 lighting = vec3(-1,5,-3);
    vec3 lightPert = lighting;
    int rep = 0;
    for (int i =0; i < rep;i++){
        lightPert.x += lightPert.x+1*rand(p.xy*time)/30;
        lightPert.y += lightPert.y+1*rand(p.xy*time+rand(p.xy*time))/30;
        lightPert.z += lightPert.z+1*rand(p.xy*time+rand(vec2(p.x)*time))/30;
    }
    vec3 diffusing = normalize((lightPert-norm));
    float sdc = rayShad(p,diffusing);    
    float d = sdc;
    float diffInt = 1*max(0*0.01*rand(time+vec2(p.x,p.y)),1*dot(norm,diffusing));    
    return 1*diffInt*d;
}

void main(void)
{
    bool initialize= true;
    vec2 uv = gl_FragCoord.xy/resolution*1 - .5;     
    vec3 prevFr = texture(prevFrame,uv+.5).rgb;
    float off = 0;
    int reps = 4;
    fragColor = vec4(vec3(0),1);
    if(initialize){
        vec3 cam = vec3(0.1*sin(time - 12*PI),0.1*cos(1*time)+.05,0.1);
        vec3 cam2 = vec3(0.1*sin((time-.5) - 12*PI),0.1*cos(1*(time-.5))+.01,0.1);
        vec3 dcam = cam - cam2;
        uv -= off;
        vec3 rd = vec3(uv,0.1);
        
        float R = 20*length(uv);
        float theta = atan(uv.y,uv.x);
        float fill = mod(log(R)*(cos(time)+1) + (theta)/(2*3.14)*2+time,.4);
        fill = 1*fill;
        vec3 back = vec3((1)*fill);    
        
        // Ray march
        sdCol tempSph= ray(cam, rd);
        vec3 shad = cam + rd*tempSph.sd;
        float lighting = getLighting(shad);
        
        // blurs, if you want it by adding scaled previous cam perspective
        vec2 offset = uv + .5 + off; 
        vec4 prevFr = texture(prevFrame, offset+dcam.xy);
        
        
        vec4 scene;
        if (tempSph.sd != 0){
            scene = vec4(1*(lighting*tempSph.col),1);
        }else{
            scene =vec4(0);
        }        
        if(tempSph.sd == 0){
            fragColor = .1*vec4(back,1)+.0*rand(uv + time);
        }else{
            fragColor = scene;
        }
        //fragColor += prevFr*.1;
    }if(true){//int(time)%reps == 0 || int(time)%reps == 1){
        vec2 duv = vec2(.01,.01);//*sin(time); 
        vec2 offset = uv + .5 + off;
        vec4 prevFr = texture(prevFrame, offset); 
        vec3 diff = vec3(0);
        vec3 diff2 = vec3(0);

        // numerical derivative on X
        diff.x = texture(prevFrame, offset+(duv.x)).x;
        diff.x -= texture(prevFrame, offset-(duv.x)).x;

        // numerical derivative on Y
        diff.y = texture(prevFrame, offset+(duv.y)).y;
        diff.y -= texture(prevFrame, offset-(duv.y)).y;        
        
        // convert to polar coords
        float w = 1*atan(uv.y,uv.x);
        float R = 1*length(uv);
        float dw = (diff.y*uv.x - diff.x*uv.y)/pow(R+.1,2);
        float dR = (diff.y*uv.y + diff.x*uv.x)/(R+.1);
        
        vec4 newFr = prevFr;
        // delta time
        float h = 0.005;
  
        // Define how each color changes with polar coordinate
        float drdw = (1)*(.5*prevFr.x+-.5*prevFr.y-0.01*prevFr.z);
        float drdR = (1)*(-.4*prevFr.x+-.3*prevFr.y-0.02*prevFr.z);
        
        float dgdR = (1)*(.5*prevFr.x+.5*prevFr.y+.002*prevFr.z);
        float dgdw = (1)*(-.9*prevFr.x+.9*prevFr.y-.003*prevFr.z);
        
        float dbdw = (1)*(.4*prevFr.x+-.01*prevFr.y+0.01*prevFr.z);
        float dbdR = (1)*(.2*prevFr.x-.1*prevFr.y+.003*prevFr.z);

        // Make equilibrium points over angle and magnitude for each color 
        float hdwdrw = h*dw*(drdw - (w + PI/2)/4);
        float hdRdrR = h*dR*(drdR - (R + 3)/4);
        
        float hdwdgw = h*dw*(dgdw - (w + PI/6)/4);
        float hdRdgR = h*dR*(dgdR - (R + 2)/4);

        float hdwdbw = h*dw*(dbdw - (w + 2*PI)/4);
        float hdRdbR = h*dR*(dbdR - (R + 1)/4);

        // integrate
        newFr.x += hdwdrw;
        newFr.x += hdRdrR;
        
        newFr.y += hdwdgw;
        newFr.y += hdRdgR;
        
        newFr.z += hdwdbw;
        newFr.z += hdRdbR;
                       
        bool diffuse = false;
        bool new = false;
        if(diffuse && !new){
            fragColor /= 10;
            fragColor += (newFr*.95);//vec4(newFr,1);//.25*new;
        }
        if(diffuse && new){
            fragColor = (newFr);//vec4(newFr,1);//.25*new;
        }        
    }
}
