import java.util.Random;
// TRIANGLE MESH
class MESH {
  // VERTICES
  int nv=0, maxnv = 1000;
  pt[] G = new pt [maxnv];
  // TRIANGLES
  int nt = 0, maxnt = maxnv*2;
  boolean[] isInterior = new boolean[maxnv];
  // CORNERS
  int c=0;    // current corner
  int nc = 0;
  int[] V = new int [3*maxnt];
  int[] O = new int [3*maxnt];
  // current corner that can be edited with keys
  MESH() {for (int i=0; i<maxnv; i++) G[i]=new pt();rand = new Random(25);};
  void reset() {nv=0; nt=0; nc=0;}                                                  // removes all vertices and triangles
  void loadVertices(pt[] P, int n) {nv=0; for (int i=0; i<n; i++) addVertex(P[i]);}
  void writeVerticesTo(pts P) {for (int i=0; i<nv; i++) P.G[i].setTo(G[i]);}
  void addVertex(pt P) { G[nv++].setTo(P); }                                             // adds a vertex to vertex table G
  void addTriangle(int i, int j, int k) {V[nc++]=i; V[nc++]=j; V[nc++]=k; nt=nc/3; }     // adds triangle (i,j,k) to V table

  // CORNER OPERATORS
  int t (int c) {int r=int(c/3); return(r);}                   // triangle of corner c
  int n (int c) {int r=3*int(c/3)+(c+1)%3; return(r);}         // next corner
  int p (int c) {int r=3*int(c/3)+(c+2)%3; return(r);}         // previous corner
  pt g (int c) {return G[V[c]];}                             // shortcut to get the point where the vertex v(c) of corner c is located

  boolean nb(int c) {return(O[c]!=c);};  // not a border corner
  boolean bord(int c) {return(O[c]==c);};  // not a border corner

  pt cg(int c) {return P(0.6,g(c),0.2,g(p(c)),0.2,g(n(c)));}   // computes offset location of point at corner c

  int ncars=0, maxncars = 1000;
  int[] cars = new int[maxncars];
  float[] locals = new float[maxncars];
  boolean[] direct = new boolean[maxncars];
  pt[] depart = new pt[maxncars];
  pt[] arrive = new pt[maxncars];
  boolean success[] = new boolean[maxncars];
  Random rand;
  //int[] road = new int[maxnv];
  float carDist = 0.3;

  int nroads = 0;
  int[] begin = new int[maxnv*10];
  int[] end = new int[maxnv*10];
  int[] roadDirects = new int[maxnv*10];
  //int[] findRoads = new int[maxnv*6];
  int[] transport = new int[maxnv*6];

  void registerRoad(){
    nroads = 0;
    for (int i = 0; i < maxnv; i++){
      transport[i] = 0;
    }
    for (int i = 0; i < nc-3; i+=3){
      for (int j = i+3; j < nc; j+=3){
        int match = 0;
        for(int k = 0; k < 3; k ++){
          if(V[i+k]==V[j] || V[i+k]==V[j+1] || V[i+k]==V[j+2]) match += 1;
          //if(V[i+k]==V[j]) match += 1;
        }
        if(match == 2){
          begin[nroads] = i;
          end[nroads] = j;
          roadDirects[nroads] = 0;
          nroads ++;
        }
      }
    }
  }

  void drawRoads(){
    fill(red);
    for(int i = 0; i < nroads; i++){
      if (roadDirects[i] == 1) arrow(triCircumcenter(begin[i]), triCircumcenter(end[i]), 10);
      if (roadDirects[i] == 2) arrow(triCircumcenter(end[i]), triCircumcenter(begin[i]), 10);
    }
  }

  // book road
  void currentRoad(int car, boolean dir){
    int end1, end2;
    end1 = car;
    if (dir) end2 = u(car);
    else end2 = s(car);
    for (int i = 0; i < nroads; i++){
      if(getRoad(end1, end2, i)==1) {
        roadDirects[i] = 1;
        transport[t(end1)]++;
        transport[t(end2)]--;
      }
      if(getRoad(end1, end2, i)==2) {
        roadDirects[i] = 2;
        transport[t(end1)]++;
        transport[t(end2)]--;
      }
    }
  }

  void nextRoad(int car, boolean dir){
    int end2, end3=-1, end4=-1;
    if (dir) {
      end2 = u(car);
      if (p(u(car)) != O[p(u(car))]) end3 = u(u(car));
      if (u(car) != O[u(car)]) end4 = s(p(u(car)));
    }
    else {
      end2 = s(car);
      if (n(s(car)) != O[n(s(car))]) end3 = s(s(car));
      if (s(car) != O[s(car)]) end4 = u(n(s(car)));
    }

    for(int i = 0; i < nroads; i ++){
      if(end4 != -1 && getRoad(end2, end4, i)==1 && roadDirects[i] != 2) {
        roadDirects[i] = 1;
        transport[t(end2)]++;
        transport[t(end4)]--;
      }
      if(end4 != -1 && getRoad(end2, end4, i)==2 && roadDirects[i] != 1) {
        roadDirects[i] = 2;
        transport[t(end2)]++;
        transport[t(end4)]--;
      }
    }

    for(int i = 0; i < nroads; i ++){
      if(end3 != -1 && getRoad(end2, end3, i)==1 && roadDirects[i] != 2) {
        roadDirects[i] = 1;
        transport[t(end2)]++;
        transport[t(end3)]--;
      }
      if(end3 != -1 && getRoad(end2, end3, i)==2 && roadDirects[i] != 1) {
        roadDirects[i] = 2;
        transport[t(end2)]++;
        transport[t(end3)]--;
      }
    }
  }

  int getRoad(int end1, int end2, int i){
    if(t(end1) == t(begin[i]) && t(end2) == t(end[i])) return 1;//{print("find", "\n"); roadDirects[nroads] = 1;arrow(triCircumcenter(begin[i]), triCircumcenter(end[i]), 20);}
    if(t(end2) == t(begin[i]) && t(end1) == t(end[i])) return 2;//{print("find", "\n"); roadDirects[nroads] = 2;arrow(triCircumcenter(end[i]), triCircumcenter(begin[i]), 20);}
    return 0;
  }

  void drawTrans(){
    fill(pink);
    for (int i = 0; i < nc; i+=3){
      if (transport[t(i)] < 0) show(triCircumcenter(i), max(abs(25*transport[t(i)]),75));
    }
  }

  // detect conflict
  boolean check(){
    for (int i = 0; i < nc; i += 3){
      if (transport[t(i)] < 0){
        int count = 0;
        for (int r = 0; r < nroads; r++){
          if (t(begin[r]) == t(i) || t(end[r]) == t(i)) count ++;
        }
        if (count == 2){
          return false;
        }
        if (count == 3 && transport[t(i)] == -2){
          return false;
        }
      }
    }
    return true;
  }

  // solve conflict
  void roadCondition(){
    for (int i = 0; i < nc; i += 3){
      if (transport[t(i)] < 0){
        int count = 0;
        for (int r = 0; r < nroads; r++){
          if (t(begin[r]) == t(i) || t(end[r]) == t(i)) count ++;
        }
        if (count == 2){
          for (int r = 0; r < nroads; r++){
            if (t(begin[r]) == t(i) && roadDirects[r] == 0){
              roadDirects[r] = 1;
              transport[t(begin[r])]++;
              transport[t(end[r])]--;
            }
            if (t(end[r]) == t(i) && roadDirects[r] == 0){
              roadDirects[r] = 2;
              transport[t(end[r])]++;
              transport[t(begin[r])]--;
            }
          }
        }
        if (count == 3 && transport[t(i)] == -2){
          for (int r = 0; r < nroads; r++){
            if (t(begin[r]) == t(i) && roadDirects[r] == 0){
              roadDirects[r] = 1;
              transport[t(begin[r])]++;
              transport[t(end[r])]--;
            }
            if (t(end[r]) == t(i) && roadDirects[r] == 0){
              roadDirects[r] = 2;
              transport[t(end[r])]++;
              transport[t(begin[r])]--;
            }
          }
        }
        if (count!=3 &&count!=2) print("error","\n");
      }
    }
  }

  boolean arriveCheck(int idx){
    pt A = findCar(cars[idx], locals[idx], direct[idx]);
    pt B = arrive[idx];
    if (norm(V(A,B)) < 100) return true;
    return false;
  }

  void updateCars(){
    float speed = 5;
    pt A, B, C, Car;
    for (int j = 0; j < nroads; j ++) roadDirects[j] = 0;


    for (int i = 0; i < ncars; i++){
      if (success[i]) continue;
      currentRoad(cars[i], direct[i]);
      if (!check()) roadCondition();
      if (!check()) roadCondition();
    }

    for (int i = 0; i < ncars; i++){
      if (success[i]) continue;
      nextRoad(cars[i], direct[i]);
      if (!check()) roadCondition();
      if (!check()) roadCondition();
    }

    for (int i = 0; i < ncars; i++){
      if (success[i]) continue;

      B = triCircumcenter(cars[i]);
      if (direct[i]){
        C = P(triCircumcenter(cars[i]), 0.5, triCircumcenter(u(cars[i])));
        A = P(triCircumcenter(cars[i]), 0.5, triCircumcenter(s(cars[i])));
      }else{
        A = P(triCircumcenter(cars[i]), 0.5, triCircumcenter(u(cars[i])));
        C = P(triCircumcenter(cars[i]), 0.5, triCircumcenter(s(cars[i])));
      }
      Car = B3(A,B,C,locals[i]);
      pt tail = B3(A,B,C,locals[i]-0.05);
      //fill(lime);
      //pillar(Car, 300, 20);
      fill(arrive[i].x % 255, arrive[i].y % 255, arrive[i].x % 255);
      show(depart[i], 15);
      arrow(arrive[i], P(arrive[i],V(0,0,100)), 10);
      hpillar(Car, tail, 40, 25, 25);

      float tmp = locals[i] + speedControl(speed, locals[i], A, B, C);
      if (tmp > 1){

        // if(direct[i]){
        //   if (rand.nextInt(2)==0 && p(u(cars[i])) != O[p(u(cars[i]))] && !conflictWhenTurn(u(cars[i]), u(u(cars[i])))){
        //     cars[i] = u(cars[i]);
        //     locals[i] = tmp - 1;
        //   }else
        //   {
        //     if(u(cars[i]) != O[u(cars[i])] && !conflictWhenTurn(p(u(cars[i])), s(p(u(cars[i]))))){
        //       cars[i] = p(u(cars[i]));
        //       direct[i] = !direct[i];
        //       locals[i] = tmp - 1;
        //     }else if(p(u(cars[i])) != O[p(u(cars[i]))] && !conflictWhenTurn(u(cars[i]), u(u(cars[i])))) {
        //       cars[i] = u(cars[i]);
        //       locals[i] = tmp - 1;
        //     }
        //   }
        //
        // }else{
        //   if (rand.nextInt(2)==0 && n(s(cars[i])) != O[n(s(cars[i]))] && !conflictWhenTurn(s(cars[i]), s(s(cars[i])))){
        //     cars[i] = s(cars[i]);
        //     locals[i] = tmp - 1;
        //   }else
        //   {
        //     if(s(cars[i]) != O[s(cars[i])] && !conflictWhenTurn(n(s(cars[i])), u(n(s(cars[i]))))){
        //       cars[i] = n(s(cars[i]));
        //       direct[i] = !direct[i];
        //       locals[i] = tmp - 1;
        //     }else if(n(s(cars[i])) != O[n(s(cars[i]))] && !conflictWhenTurn(s(cars[i]), s(s(cars[i])))) {
        //       cars[i] = s(cars[i]);
        //       locals[i] = tmp - 1;
        //     }
        //   }
        // }

        int left, right, togo;
        left = nextLeft(cars[i], direct[i]);
        right = nextRight(cars[i], direct[i]);
        togo = toGo(left, right, arrive[i]);

        if (togo == -1){
          cars[i] = left;
          direct[i] = true;
          locals[i] = tmp - 1;
        }else if (togo == 1){
          cars[i] = right;
          direct[i] = false;
          locals[i] = tmp - 1;
        }
        currentRoad(cars[i], direct[i]);
        if (!check()) roadCondition();
        if (!check()) roadCondition();
      }
      else locals[i] = tmp;

      success[i] = arriveCheck(i);

    }
  }

  int nextLeft(int car, boolean direct){
    if (direct){
      if (p(u(car)) == O[p(u(car))] || conflictWhenTurn(u(car), u(u(car))) ) return -1;
      else return u(car);
    }else{
      if (s(car) == O[s(car)] || conflictWhenTurn(n(s(car)), u(n(s(car))))) return -1;
      else return n(s(car));
    }
  }

  int nextRight(int car, boolean direct){
    if (direct){
      if(u(car) == O[u(car)] || conflictWhenTurn(p(u(car)), s(p(u(car))))) return -1;
      else return p(u(car));
    }else{
      if(n(s(car)) == O[n(s(car))] || conflictWhenTurn(s(car), s(s(car)))) return -1;
      else return s(car);
    }
  }

  int toGo(int left, int right, pt arrive){
    if (left == -1 && right == -1){
      return 0;
    }else if (left == -1){
      return 1;
    }else if (right == -1){
      return -1;
    }else{
      float tmp1, tmp2;
      tmp1 = norm(V(findCar(left, 0.5, true), arrive));
      tmp2 = norm(V(findCar(right, 0.5, true), arrive));
      if (tmp1 > tmp2) return 1;
      else return -1;
    }
  }

  void addCar(){
    cars[ncars] = rand.nextInt(nc);
    while(!isInterior[V[cars[ncars]]]){
      cars[ncars] = rand.nextInt(nc);
    }
    locals[ncars] = 1. * rand.nextInt(100) / 100;
    direct[ncars] = (rand.nextInt(2)==0);

    while(conflictWhenAdd(cars[ncars], locals[ncars], direct[ncars])){
      cars[ncars] = rand.nextInt(nc);
      while(!isInterior[V[cars[ncars]]]){
        cars[ncars] = rand.nextInt(nc);
      }
      locals[ncars] = 1. * rand.nextInt(100) / 100;
      direct[ncars] = (rand.nextInt(2)==0);
    }

    depart[ncars] = findCar(cars[ncars],locals[ncars], direct[ncars]);

    int car = rand.nextInt(nc);
    while(!isInterior[V[car]]){
      car = rand.nextInt(nc);
    }
    float local = 1. * rand.nextInt(100) / 100;
    arrive[ncars] = findCar(car,local, true);

    success[ncars] = false;

    ncars ++;
  }

  pt findCar(int car, float local, boolean dir){
    pt A, B, C;
    B = triCircumcenter(car);
    if (dir){
      C = P(triCircumcenter(car), 0.5, triCircumcenter(u(car)));
      A = P(triCircumcenter(car), 0.5, triCircumcenter(s(car)));
    }else{
      A = P(triCircumcenter(car), 0.5, triCircumcenter(u(car)));
      C = P(triCircumcenter(car), 0.5, triCircumcenter(s(car)));
    }
    return B3(A,B,C,local);
  }

  float speedControl(float speed, float t, pt A, pt B, pt C){
    float vx = 2*t*C.x - (4*t-2)*B.x + (2*t-2)*A.x;
    float vy = 2*t*C.y - (4*t-2)*B.y + (2*t-2)*A.y;
    return speed / (float)Math.sqrt(vx*vx + vy*vy);
  }

  boolean conflictWhenAdd(int car, float local, boolean dir){
    int end1, end2;
    end1 = car;
    if (dir) end2 = u(car);
    else end2 = s(car);
    for (int i = 0; i < nroads; i++){
      if(getRoad(end1, end2, i)==1 && roadDirects[i] ==2) {
        return true;
      }else if (getRoad(end1, end2, i)==2 && roadDirects[i] ==1){
        return true;
      }
    }

    for (int i = 0; i < ncars; i++){
      if (success[i]) continue;
      if(cars[i] == car){
        if(direct[i] == dir){
          if(abs(local - locals[i])<carDist) return true;
        }else return true;
      }
    }
    return false;
  }

  boolean conflictWhenTurn(int end2, int end4){
    for(int i = 0; i < nroads; i ++){
      if(getRoad(end2, end4, i)==1 && roadDirects[i] != 2) return false;
      if(getRoad(end2, end4, i)==2 && roadDirects[i] != 1) return false;
    }
    return true;
  }




  // CORNER ACTIONS CURRENT CORNER c
  void next() {c=n(c);}
  void previous() {c=p(c);}
  void opposite() {c=o(c);}
  void left() {c=l(c);}
  void right() {c=r(c);}
  void swing() {c=s(c);}
  void unswing() {c=u(c);}
  void printCorner() {println("c = "+c);}



  // DISPLAY
  void showCurrentCorner(float r) { if(bord(c)) fill(red); else fill(dgreen); show(cg(c),r); };   // renders corner c as small ball
  void showDirect(){
    //for (int c = 0; c < nc; c++){
    //  if (isInterior[V[c]]){
    //    if (road[c] == 0) {fill(white); show(P(0.6,g(c),0.2,g(p(c)),0.2,g(n(c))), 10);}
    //    if (road[c] == 1) {fill(orange); show(P(0.62,g(c),0.2,g(p(c)),0.2,g(n(c))), 10);}
    //    if (road[c] == 2) {fill(black); show(P(0.58,g(c),0.2,g(p(c)),0.2,g(n(c))), 10);}
    //  }
    //}
  }
  void showEdge(int c) {beam( g(p(c)),g(n(c)),rt ); };  // draws edge of t(c) opposite to corner c
  void showVertices(float r) // shows all vertices green inside, red outside
    {
    for (int v=0; v<nv; v++)
      {
      if(isInterior[v]) fill(green); else fill(red);
      show(G[v],r);
      }
    }
  void showInteriorVertices(float r) {for (int v=0; v<nv; v++) if(isInterior[v]) show(G[v],r); }                          // shows all vertices as dots
  void showTriangles() { for (int c=0; c<nc; c+=3) show(g(c), g(c+1), g(c+2)); }         // draws all triangles (edges, or filled)
  void showEdges() {for (int i=0; i<nc; i++) showEdge(i); };         // draws all edges of mesh twice

  void triangulate()      // performs Delaunay triangulation using a quartic algorithm
   {
     c=0;                   // to reset current corner
     // **01 implement it
     int i, j, k, m;
     for (i = 0; i < nv; i++){
       for (j = i+1; j < nv; j++){
         for (k = j+1; k < nv; k++){
           boolean good = true;
           pt O = CircumCenter(G[i],G[j],G[k]);
           for (m = 0; m < nv; m++){
             if (m == i || m == j || m == k) continue;
             if (norm(V(O, G[k])) > norm(V(O, G[m]))) {
               good = false;
               break;
             }
           }
           if (good) {
             vec Normal = cross(V(G[i],G[j]), V(G[j],G[k]));
             //print(triThickness(G[i], G[j], G[k]));
             if (!isFlatterThan(G[i], G[j], G[k], 10)){
                 if (Normal.z < 0) addTriangle(i, k, j);
                 else addTriangle(i, j, k);
             }
           }
         }
       }
     }
   }

   boolean isFlatterThan( pt A, pt B, pt C, float r) {return triThickness(A, B, C) < r;}


  void computeO() // **02 implement it
    {
    // **02 implement it
      int i, j;
      for (i = 0; i < nc; i++) O[i] = i;
      for (i = 0; i < nc; i++){
        for (j = i+1; j < nc; j++){
          //print(i, V[i], j, V[j], V[n(i)], V[p(j)], V[n(j)], V[p(i)], '\n');
          if ((V[n(i)] == V[p(j)]) && (V[n(j)] == V[p(i)])){
            O[i] = j; O[j] = i;
          }
        }
      }
    }

  void showBorderEdges()  // draws all border edges of mesh
    {
    // **02 implement;
      for (int i = 0; i < nc; i++){
        if (O[i] == i) cylinderSection(g(p(i)), g(n(i)), 20);
      }
    }

  void showNonBorderEdges() // draws all non-border edges of mesh
    {
    // **02 implement
      for (int i = 0; i < nc; i++){
        if (O[i] != i) beam(g(p(i)),g(n(i)),10);
      }
    }

  int countBorders(){
    int count = 0;
    for (int i = 0; i < nc; i++){
      if (O[i]==i) count++;
    }
    return count;
  }

  void classifyVertices()
    {
    // **03 implement it
      for(int i = 0; i < nv; i++) isInterior[i] = true;
      for(int i = 0; i < nc; i++){
        if (bord(i)){
          isInterior[V[p(i)]] = false;
          isInterior[V[n(i)]] = false;
        }
      }
    }

  void smoothenInterior(float t) { // even interior vertiex locations
      pt[] Gn = new pt[nv];
      // **04 implement it
      for (int v = 0; v < nv; v++){
        if (!isInterior[v]) continue;
        Gn[v] = P(0,0,0);
        float count = 0;
        for(int i = 0; i < nc; i++){
          if (V[i] == v){
            float weight = norm(V(g(i),g(n(i))));
            //float weight = 0;
            //for(int j = 0; j < nc; j++){
            //  if (V[j] == V[n(i)]){
            //    weight += norm(V(g(n(j)),g(n(i))));
            //  }
            //}
            Gn[v] = P(count/(count+weight), Gn[v], weight/(count+weight), g(n(i)));
            count += weight;
          }
        }
        //print(v, count, '\n');
      }
      for (int v=0; v<nv; v++) if(isInterior[v]) G[v].translateTowards(t,Gn[v]);
    }


   // **05 implement corner operators in Mesh
  int v (int c) {return V[c];}                                // vertex of c
  int o (int c) {return O[c];}                                // opposite corner
  int l (int c) {return O[n(c)];}                             // left
  int s (int c) {return n(l(c));}                             // left
  int u (int c) {return p(r(c));}                             // left
  int r (int c) {return O[p(c)];}                             // right

  void showOpposites(){
    for(int i = 0; i < nc; i++){
      if (O[i] > i) {
        //print(i, O[i], '\n');
        pt A = g(i);
        pt C = g(O[i]);
        pt B = P(g(n(i)), 0.5, g(p(i)));
        for (float j = 0; j < 1; j+=0.1) beam(N3(A,B,C,j), N3(A,B,C,j+0.1), 0.1);
      }
    }
  }

  pt N3(pt A, pt B, pt C, float t){
    pt P = P(P(A, 2*t, B), t, P(B, 2*t-1, C));
    return P;
  }

  void showVoronoiEdges() // draws Voronoi edges on the boundary of Voroni cells of interior vertices
    {
    // **06 implement it
      for (int i=0; i<nc; i++){
        if(O[i] > i){
          beam(triCircumcenter(i), triCircumcenter(O[i]), 3);
        }
      }
    }

  void showArcs() // draws arcs of quadratic B-spline of Voronoi boundary loops of interior vertices
    {
    // **06 implement it
      for (int i=0; i<nc; i++){
        if (V[s(i)] == V[i] && V[u(i)] == V[i]){
          pt B = triCircumcenter(i);
          pt A = triCircumcenter(s(i)); A = P(B, 0.5, A);
          pt C = triCircumcenter(u(i)); C = P(B, 0.5, C);
          for (float j = 0; j < 1; j+=0.1) beam(B3(A,B,C,j), B3(A,B,C,j+0.1), 1);
        }
      }
    }               // draws arcs in triangles

  pt B3(pt A, pt B, pt C, float t){
    pt P = P(P(A,t,B), t, P(B,t,C));
    return P;
  }

  void drawVoronoiFaceOfInteriorVertices(){
    float dc = 1./(nv - 1);
    for (int v = 0; v < nv; v++)
    {
      if (isInterior[v]) {
        fill(dc*255*v, dc*255*(nv-v), 200);
        //fill(black);stroke(black);
        drawVoronoiFaceOfInteriorVertices(v);
      }
    }
  }

  void drawVoronoiFaceOfInteriorVertices(int v){
    int center = 0;
    for (int i = 0; i < nc; i++){
      if(V[i] == v) {center = i; break;}
    }
    int neighbor = center;

    while (true) {
      beginShape(TRIANGLE_FAN);
        vert(g(center));
        vert(triCircumcenter(neighbor));
        vert(triCircumcenter(s(neighbor)));
      endShape();
      //print("111/n");
      neighbor = s(neighbor);
      if (neighbor == center) break;
    }

  }

  pt triCenter(int c) {return P(g(c),g(n(c)),g(p(c))); }  // returns center of mass of triangle of corner c
  pt triCircumcenter(int c) {return CircumCenter(g(c),g(n(c)),g(p(c))); }  // returns circumcenter of triangle of corner c


  } // end of MESH
