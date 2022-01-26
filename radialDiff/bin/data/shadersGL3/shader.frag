#version 150
#ifdef GL_ES
    precision mediump float;
#endif
uniform float time;
uniform float h;
uniform vec2 resolution;
uniform vec2 mouse;
uniform int diffuse;
uniform int stop;
uniform vec3 heat;
uniform vec3 camPos;
uniform vec3 mcamPos;
uniform vec3 scamPos;

uniform sampler2D texture0;
uniform sampler2D texture1;
uniform sampler2D texture2;
uniform sampler2D texture3;
uniform sampler2D video;
uniform sampler2D prevFrame;
uniform sampler2D prevPass;

#define PI 3.141592653
in VertexData
{
    vec4 v_position;
    vec3 v_normal;
    vec2 v_texcoord;
} inData;

struct sdCol{
    float sd;
    vec3 col;
};

vec3 matMult(mat3 Map, vec3 vec){
    vec3 temp = vec;
    temp.x = dot(Map[0],temp);
    temp.y = dot(Map[1],temp);
    temp.z = dot(Map[2],temp);
    return temp;
}

float rand(vec2 co){
    return fract(sin(dot(co, vec2(12.9898, 78.233))) * 43758.5453);
}

out vec4 fragColor;
float RP = exp(cos(time)-1);
float thetaP = mod(time,2*3.14);
float phiP = mod(2*time-2*3.14,2*3.14);
vec3 pos = vec3(RP*sin(thetaP)*cos(phiP),RP*sin(thetaP)*sin(phiP)+.3,1.-.1*sin(time));
vec3 pos2 = vec3(1*cos(time),1*sin(cos(sin(time)))+.1,1);

sdCol sdGround(vec3 p){
    sdCol sdG;
    sdG.sd = 0;
    mat3 view = mat3(.1*cos(.2*time),-.1*sin(.2*time),0,//sin(time),
                .1*sin(.2*time),.1*cos(.2*time),0,
                0,0,1);
    // Mirror over x and y dimensions
    // Threshold on Z dimension
    vec3 sdd = vec3(cos(((p.x))),(p.y),(((-.01*p.z))));
    sdG.sd = (sdd.y)*0.5+(sdd.x + sdd.z) +  .5*sin(2.0 * p.x)*cos(3*p.x) + 0*cos(p.z);//(atan(p.y-1,p.x));    
    //sdG.sd += sin(5*p.x)*sin(5*p.y)-sin(p.z);
    sdG.col = vec3(1*sin(p.x)*+1.1,(p.y)+1.1,exp(-p.z));//sin(time)*cos(p.x),0.01*cos(p.y),exp(-20*p.y));//3*sin(p.x)*cos(time),0.1*cos(p.y)*sin(3*time),cos(sin(p.z))*cos(time-100));
    return sdG;
}

sdCol sdSpheres(vec3 p, vec3 cent, float rad, vec3 col){
    sdCol sdS;
    sdS.sd = length(p -cent) - rad;
    sdS.col = col;
    return sdS;
}

sdCol sdRectPrism(vec3 p, vec3 cent, vec3 dim, vec3 col){
    sdCol sdS;
    //sdS.sd = length(p - cent);
    dim = dim;
    vec3 surfDist = abs((mat3(1,0, -1,0,1,0,1,0,1))*p - cent) - (1)*dim;
    sdS.sd = length(max(surfDist,0)) + min(max(surfDist.x, max(surfDist.y,surfDist.z)),0);
    
    sdS.col = col;
    return sdS;
}

sdCol sdTor(vec3 p, vec3 cent, float R1, float R2, vec3 col){
//    vec3 q = p - clamp(p,vec3(p,p.x-1,0,p.z-1),vec3(p.x+1,0,pz+1));
    sdCol tor;
    vec3 center = p - cent;
    vec2 q = vec2(length(center.xy) - R1, center.z);
    tor.sd = length(q) - R2;
    tor.col = col;
    return tor;
}

sdCol minWithCol(sdCol obj1, sdCol obj2){
    if (obj1.sd < obj2.sd){
        return obj1;
    }else{
        return obj2;
    }
}

sdCol addShad(sdCol obj1,sdCol obj2){
    sdCol temp = obj1;
    temp.sd = obj1.sd*obj2.sd;
    temp.col = obj1.col/(obj2.col);
    return temp;
}

sdCol mapofWorld(vec3 p){

    sdCol ground = sdGround(p);
    float Rp = length(p);
    float theta = atan(p.y,p.x);
    float phi = atan(length(p.xy),p.z);
    float rad1 = .3;// + Rp*sin(5*theta)*sin(5.0*phi)*0.1;
    float rad2 = .2;// + Rp*sin(5*theta)*sin(5.0*phi)*0.25;;// +sin(phi*2)*sin(2*theta)*.25;//+.1*(cos(3*phi*theta*Rp)+.2); 
    vec3 colsph1 = 4*vec3(1*p.z,2*p.y,0.1*cos(p.x));
    sdCol Sphere1 = sdSpheres(p, pos, rad1, colsph1);    
    sdCol Sphere2 = sdSpheres(p, pos2, rad2, 3*vec3(cos(time)*cos(sin(p.z))*.1,.4*sin(p.z),cos(0)));
    sdCol tor = sdTor(p, pos, .4*rad1, .1*rad1, colsph1);
    tor.sd += 0;
    sdCol Rect = sdRectPrism(p, pos2, vec3(.1), 3*vec3(0*cos(time)*cos(sin(p.z))*.1,0*.4*sin(p.z),cos(time)));    
    sdCol tempCol = minWithCol(tor,Sphere2);
    //tempCol = minWithCol(tempCol, Rect);
    /*sdCol ground;
    ground.sd = p.y/1+1;
    ground.col = vec3(1,0,0);*/
    tempCol = minWithCol(tempCol, ground);
    return tempCol;
}


const float MAXIMUM_TRACE_DISTANCE = 100.0;

vec3 calcNorm(vec3 p){
    const vec3 small_step = vec3(0.001, 0.00, 0.00);
    float gradient_x = mapofWorld(p + small_step.xyy).sd - mapofWorld(p - small_step.xyy).sd;
    float gradient_y = mapofWorld(p + small_step.yxy).sd - mapofWorld(p - small_step.yxy).sd;
    float gradient_z = mapofWorld(p + small_step.yyx).sd - mapofWorld(p - small_step.yyx).sd;

    vec3 normal = vec3(gradient_x, gradient_y, gradient_z);

    return normalize(normal);
}

sdCol ray(vec3 r, vec3 rd){
    vec3 ro = r;  
    float dO = 0;
    float tol =1e-1;    
    sdCol dist;
    dist.sd = -1;
    dist.col = vec3(0);
    float stepUp = 0;
    for (int i = 0; i < 32; i++){
        ro += ro + rd*(dO);
        sdCol dist = mapofWorld(ro);
        dO += dist.sd;
        if (dist.sd - stepUp < 0.001){
            dO += .05;
        }
        if (dist.sd < tol){
            dist.sd = dO;
            return dist;//ro +dist.col;
        }
        if (dO >= 1e10){
            dist.sd = -1;
            break;
        }
    }

    dist.sd = -1;
    return dist;//back/10;
}

float rayShad(vec3 r, vec3 rd){
    vec3 ro = r;  
    float dO = 0;
    float tol = 1e-10;
    float tmint = 0;
    float tmaxt = 10;
    float res = 1;
    int k = 128;
    float ph = 1e10;
    for (float t = tmint; t < tmaxt;){
        ro += ro + rd*(t);
        sdCol dist = mapofWorld(ro);
        float y = dist.sd*dist.sd/(2*ph);
        float d = sqrt(dist.sd*dist.sd - y*y);
        t += dist.sd;
        res = min(res,k*dist.sd/max(0,t-y));
        ph = dist.sd;
        if (dist.sd < tol){            
            return 0;
        }
    }
    return res;
}

float getLighting(vec3 p, vec2 uv){
    vec3 norm = calcNorm(p);
    int lights = 2;
    vec3 lighting[2] = vec3[2](vec3(-1,1,-1),vec3(10,10,1));//-1*sin(time)*sin(time),1*cos(time)*cos(time),-8*sin(time)*sin(time));//1*sin(time+PI*3)*.1+2,1*sin(time)+1,-sin(time)-.5);//1*sin(time+PI*3)*5-2.5,1*sin(time)-1,-sin(time)-1);
    float diffInt =0;
    for (int i = 0; i < lights; i++){
        vec3 lightPert = lighting[i];
        vec3 diffusing = 1*normalize((lightPert-norm));
        vec3 rr = normalize(cross(vec3(0,1,0), diffusing));
        vec3 ru = normalize(cross(diffusing, rr));
        float fov = 1;
        vec3 diff = normalize((uv.x+sin(PI/2)) * rr + (uv.y+cos(0)) * ru + diffusing*fov);
        float sdc = rayShad(p+norm,diff);    
        float d = sdc;
        diffInt += 1*(d)*max(0,dot(norm,diff));
    }
    return 1*diffInt/lights;
}

vec4 integrate(vec4 prev){
    float sigma = 10;
    float rho = 28;
    float beta = 8/3;
    vec3 diff = vec3(0);
    /*
    diff.x = -prev.y - prev.z;//sigma*(prev.y - prev.x);
    diff.y = prev.x + 0.2*prev.y;//prev.x*(rho-prev.z)-prev.y;    
    diff.z = 0.2 + prev.z*(prev.x - 5.7);//prev.x*(prev.y)-beta*prev.z;
    */
    diff.x = sigma*(prev.y - prev.x);
    diff.y = prev.x*(rho-prev.z)-prev.y;    
    diff.z = prev.x*(prev.y)-beta*prev.z;    
    float dt = 1e-3;
    float dn = 1e-4;
    vec3 o = vec3(0);
    o.x = prev.x + diff.x*dt + dn*rand(vec2(diff.z*time));
    o.y = prev.y + diff.y*dt + dn*rand(vec2(prev.x*time));
    o.z = prev.z + diff.z*dt + dn*rand(vec2(diff.x*time));    
    return vec4(o,1);
}

void main(void)
{    
    vec2 uv = 2*(gl_FragCoord.xy/resolution*1 - .5);     

    float off = (sin(time)*tanh(length(camPos.xy))*((uv.y))) + (sin(time)*sin(uv.x));
    int reps = 4;
    fragColor = vec4(vec3(0),1);
    vec3 cam = 1/scamPos*(camPos - mcamPos)+vec3(0,0,.1);//0.1*sin(time - 12*PI),0.1*cos(1*time)+.05,0.1);//0.25*sin((0*time)),0.05*sin(0*time),0.01*cos(sin(2*time)));
//    vec3 cam2 = vec3(0.1*sin((time-.5) - 12*PI),0.1*cos(1*(time-.5))+.01,0);//0.25*sin((0*time)),0.05*sin(0*time),0.01*cos(sin(2*time)));
//    vec3 dcam = cam - cam2;
    uv -= off;//sin(tanh(uv.y))*sin(time*2) + tanh(sin(time)*(uv.x));
    uv *= vec2(1,-1);
    vec3 rd = normalize(-(cam - vec3(uv,1)));//*sin(time));
    vec3 rr = normalize(cross(vec3(0,1,0), rd));
    vec3 ru = normalize(cross(rd,rr));
    float fPersp = 1;

    rd = normalize((uv.x)*rr + (uv.y)*ru + rd*fPersp);        
    float R = 1*length(uv);
    float theta = atan(uv.y,uv.x);
    float fill = exp(-5*mod(log(R)*(cos(time)+1) + (theta)/(2*3.14)*2+time,.4));
    fill = 1*fill;//exp(abs(sin(.2*time))*abs(sin(time))*fill);
    vec3 back = vec3((1)*fill);
    sdCol tempSph= ray(cam, rd);
    vec3 shad = cam + rd*tempSph.sd;
    float lighting = getLighting(shad, uv);

    // blurs
    vec2 offset = uv + 1 + off; 
//    vec4 prevFr = texture(prevFrame, offset+dcam.xy);


    vec4 scene;
    vec4 tex4 = texture(video,(uv*vec2(1,-1)+1)/2);// + 1 + off));
    if (tempSph.sd != -1){
        scene = vec4(tempSph.col*lighting,1)*1;
    }else{
        scene = vec4(back*rd.y*.6*vec3(tan(cos(time)),sin(time),cos(time)),1);        
    }
    //scene += .1*(.1+.1*sin(time))*tex4;
    
    scene.x = pow(scene.x,1/2.2);
    scene.y = pow(scene.y,1/2.2);
    scene.z = pow(scene.z,1/2.2);
    if (length(tex4.xyz) > .9){
        fragColor = scene;
    }else{
        fragColor = tex4;//200;
    }
    if(true){//int(time)%reps == 0 || int(time)%reps == 1){
        vec2 duv = vec2(.1,.1);//*sin(time); 
        uv *= vec2(1,-1);
        vec2 offset = (uv + 1 + off)/2;// + sin(tanh(uv.y))*sin(time*2) + tanh(sin(time)*(uv.x));
        vec4 prevFr = texture(prevFrame, offset); //
        vec3 diff = vec3(0);
        vec3 diff2 = vec3(0);
        diff.x = texture(prevFrame, offset+(duv.x)).x;
        diff.x -= texture(prevFrame, offset-(duv.x)).x;
        //diff.x += texture(prevFrame, offset+(duv.x)).y;
        //diff.x -= texture(prevFrame, offset-(duv.x)).y;
        diff.y = texture(prevFrame, offset+(duv.y)).y;
        diff.y -= texture(prevFrame, offset-(duv.y)).y;
        //diff.y += texture(prevFrame, offset+(duv.y)).x;
        //diff.y -= texture(prevFrame, offset-(duv.y)).x;
        
        float w = 1*atan(uv.y,uv.x);
        //float dR = length(vec2(diff.x, diff.y));
        float R = 1*length(uv);
        float Rw = mod(R + cos(w)/(2*PI),2*PI);
        float dw = (diff.y*uv.x - diff.x*uv.y)/pow(R+.5,2);//atan(diff.y,diff.x);
        float dR = (diff.y*uv.y + diff.x*uv.x)/(R+0.001);
        
        vec4 posTex = texture(texture0, offset + pos.xy);
        vec4 newFr = prevFr;
        //float h = .05;//+.0*(sin(time)+1);


        diff = (diff);
        R = length(uv);
        for (int i = 0; i < 1; i++){
        float drdw = (0)*(.5*prevFr.x+-.5*prevFr.y-0.01*prevFr.z);
        float drdR = (0)*(-.4*prevFr.x+-.3*prevFr.y-0.02*prevFr.z);
        
        float dgdR = (0)*(.5*prevFr.x+.5*prevFr.y+.002*prevFr.z);
        float dgdw = (0)*(-.9*prevFr.x+.9*prevFr.y-.003*prevFr.z);
        
        float dbdw = (0)*(.4*prevFr.x+-.01*prevFr.y+0.01*prevFr.z);
        float dbdR = (0)*(.2*prevFr.x-.1*prevFr.y+.003*prevFr.z);

        float hdwdrw = h*dw*(2*heat.x+drdw - (mod(w,2*PI) - PI/2)/4);
        float hdRdrR = h*dR*(2*heat.x+drdR -( R - 7)/4);
        
        float hdwdgw = h*dw*(2*heat.y+dgdw - (mod(w,2*PI) - PI/6)/4);
        float hdRdgR = h*dR*(2*heat.y+dgdR - (R - 3)/4);

        float hdwdbw = h*dw*(2*heat.z+dbdw - (mod(w,2*PI) - 2*PI)/4);
        float hdRdbR = h*dR*(2*heat.z+dbdR - (R - 5)/4);
                
        newFr.x += hdwdrw;
        newFr.x += hdRdrR;
        
        newFr.y += hdwdgw;
        newFr.y += hdRdgR;
        
        newFr.z += hdwdbw;
        newFr.z += hdRdbR;
        prevFr = newFr;
        }
//        bool diffuse = true;
        if(diffuse==1 && stop == 0){
            fragColor /= 10;
            fragColor += (newFr*.97);
        }
        if(diffuse==1 && stop==1){
            fragColor = (newFr);
        }
    }
}
