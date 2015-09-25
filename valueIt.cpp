#include <cstdlib>
#include <algorithm>
// #include <Rcpp.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

class Coord{
 public:
  Coord(const int x, const int y);
  Coord(const Coord & c);
  bool operator==(const Coord & rhs) const;
  int x;
  int y;
};


Coord::Coord(const int x, const int y)
  :x(x),y(y){
};

Coord::Coord(const Coord & c)
  :x(c.x),y(c.y){
}


bool Coord::operator==(const Coord & rhs) const{
  return this->x == rhs.x && this->y == rhs.y;
}

class Action{
 public:
  Action(const int x, const int y);
  int x;
  int y;
};

Action::Action(const int x, const int y)
  :x(x),y(y){
};


class Grid {
 public:
  Grid(const std::vector<double> r,
       const int n,
       const Coord d,
       const Coord g,
       const Coord s,
       const double noise,
       const int nActions,
       const std::vector<Action> actions);

  std::vector<double> r;
  int n;
  Coord d;
  Coord g;
  Coord s;
  double noise;
  int nActions;
  std::vector<Action> actions;

  int c2i(const Coord & s);

  Coord move(const Coord & s, const Action & a);
  double transProb(const Coord & s, const Action & a, const Coord & sp);
  double expReward(const Coord & s, const Action & a);
};


Grid::Grid(const std::vector<double> r,
	   const int n,
	   const Coord d,
	   const Coord g,
	   const Coord s,
	   const double noise,
	   const int nActions,
	   const std::vector<Action> actions)
  :r(r),n(n),d(d),g(g),s(s),noise(noise),nActions(nActions),actions(actions){
};



int Grid::c2i(const Coord & s){
  return s.y * d.x + s.x;
}


Coord Grid::move(const Coord & s, const Action & a){
  Coord sp(std::max(std::min(d.x-1,s.x + a.x),0),
	std::max(std::min(d.y-1,s.y + a.y),0));
  return sp;
}


double Grid::transProb(const Coord & s, const Action & a, const Coord & sp){
  if(s == g && s == sp)
    return 1.0;
  else if(s == g)
    return 0;
  else{
    // if adherence
    Coord sa = move(s,a);
    double prob = 0.0;
    if(sa == sp)
      prob += 1.0 - noise;

    // if noisy action
    int i = 0;
    for(i = 0; i < nActions; ++i){
      sa = move(s,actions[i]);
      if(sa == sp)
	prob += noise/double(nActions);
    }
    return(prob);
  }
}


double Grid::expReward(const Coord & s, const Action & a){
  if(s == g)
    return 0.0;
  else{
    // if adherence
    Coord sa = move(s,a);
    double er = r.at(c2i(sa)) * (1.0 - noise);

    // if noisy
    int i;
    for(i = 0; i < nActions; ++i){
      sa = move(s,actions[i]);
      er += r.at(c2i(sa)) * noise / double(nActions);
    }
    return(er);
  }
}

// [[Rcpp::export]]
arma::colvec
solveValueIterFast(arma::colvec v,
		   const std::vector<double> r,
		   const int dX,
		   const int dY,
		   const int gX,
		   const int gY,
		   const double noise,
		   const double sX,
		   const double sY,
		   const std::vector<int> actions,
		   const double gamma,
		   const double tol){
  int i;
  std::vector<Action> actions_;
  for(i = 0; i < int(actions.size()/2); ++i){
    actions_.push_back(Action(actions.at(2*i),
  			      actions.at(2*i+1)));
  }

  Grid g(r,dX*dY,Coord(dX,dY),Coord(gX,gY),Coord(sX,sY),
  	 noise,actions_.size(),actions_);

  std::vector<arma::colvec> eR(g.actions.size());
  std::vector<arma::mat> T(g.actions.size());

  int a,x,y,xp,yp;
  for(a = 0; a < g.nActions; ++a){
    eR.at(a).resize(g.n);
    T.at(a).resize(g.n,g.n);
    for(x = 0; x < g.d.x; ++x){
      for(y = 0; y < g.d.y; ++y){
  	Coord s(x,y);
  	eR.at(a)(g.c2i(s)) = g.expReward(s,g.actions[a]);
  	for(xp = 0; xp < g.d.x; ++xp){
  	  for(yp = 0; yp < g.d.y; ++yp){
  	    Coord sp(xp,yp);
  	    T.at(a)(g.c2i(s),g.c2i(sp)) = g.transProb(s,g.actions[a],sp);
  	  }
  	}
      }
    }
  }

  arma::colvec vp = v;
  bool cont = true;
  int s;
  while(cont){
    std::vector<arma::colvec> vA(g.nActions);
    for(a = 0; a < g.nActions; ++a){
      vA[a] = eR[a] + gamma * T[a] * v;
    }


    for(s = 0; s < g.n; ++s){
      std::vector<double> vsa(g.nActions);
      for(a = 0; a < g.nActions; ++a){
  	vsa[a] = vA[a][s];
      }
      v[s] = *std::max_element(vsa.begin(),vsa.end());
    }

    double diff = arma::norm(vp - v,2);
    cont = diff > tol;
    vp = v;
    printf("% 16.8f :: % 16.8f\n",arma::sum(v),diff);
  }
  return v;
  // return arma::colvec(1);
}
