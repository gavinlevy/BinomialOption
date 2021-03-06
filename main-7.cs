/******************************************************************************

                            Online C# Compiler.
                Code, Compile, Run and Debug C# program online.
Write your code in this editor and press "Run" button to execute it.

*******************************************************************************/

using System;
class Project3 {
  static void Main() {
    
     //bool OptionType = true; //true is call, false is put
        bool OptionKind = true; //true is american, false is european
        double Strike = 90;
        double S0 = 88; //underlying
        double ExpireTime = 10;
        double t0 = 0;
        int n = 5;
        double RiskFreeRate = 0.02;
        double ContinuousDividendYield = 0.0;
        double Volatility = 0.2;
        double Vtarget=18;
        double Vtol=0.01;

        BinomialOption p = new BinomialOption(false, OptionKind, Strike, ExpireTime, RiskFreeRate, ContinuousDividendYield, Volatility);
        BinomialOption c = new BinomialOption(true, OptionKind, Strike, ExpireTime, RiskFreeRate, ContinuousDividendYield, Volatility);
        
        //put call parity
        Console.WriteLine("Put Call Parity Check");
        Console.WriteLine("C+Ke^(-rt): " + (c.Valuation(S0, t0, n) + (Strike * (Math.Exp(-RiskFreeRate * ExpireTime)))));
        Console.WriteLine("P + S0: " + (p.Valuation(S0, t0, n) + S0));
        
        Console.WriteLine(" ");
        Console.WriteLine("The Put option value is:" + p.Valuation(S0, t0, n));
        Console.WriteLine("The Delta is:" + p.Del(S0, t0, n));
        Console.WriteLine("The Gamma is:" + p.Gam(S0, t0, n));
        Console.WriteLine("The Veg is:" + p.Veg(S0, t0, n));
        Console.WriteLine("The Theta is:" + p.Thet(S0, t0, n));
        Console.WriteLine("The Rho is:" + p.Rhoo(S0, t0, n));
        Console.WriteLine("The implied volatility is : " + p.impliedVolatility(S0, t0, Vtarget, Vtol, n));


        Console.WriteLine(" ");
        Console.WriteLine("The Call option value is:" + c.Valuation(S0, t0, n));
        Console.WriteLine("The Delta is:" + c.Del(S0, t0, n));
        Console.WriteLine("The Gamma is:" + c.Gam(S0, t0, n));
        Console.WriteLine("The Veg is:" + c.Veg(S0, t0, n));
        Console.WriteLine("The Theta is:" + c.Thet(S0, t0, n));
        Console.WriteLine("The Rho is:" + c.Rhoo(S0, t0, n));
        Console.WriteLine("The implied volatility is : " + c.impliedVolatility(S0, t0, Vtarget, Vtol, n));
  }
}

class BinomialOption
{
  private bool type;
  private bool kind;
  private double K;
  private double T;
  private double r;
  private double q;
  private double ??;

  public BinomialOption (bool OptionType, bool OptionKind, double Strike, double ExpireTime, double RiskFreeRate, double ContinuousDividendYield, double Volatility)
  {
    type = OptionType;
    kind = OptionKind;
    K = Strike;
    T = ExpireTime;
    r = RiskFreeRate;
    q = ContinuousDividendYield;
    ?? = Volatility;
  }

  public string getOptionType()
  {
      if(type == true) { return "Call Option."; }
      else             { return "Put Option."; }
  }
  public string getOptionKind()
  {
      if(kind == true) { return "American Option."; }
      else             { return "European Option."; }
  }
  public double getStrike ()                  {return K;}
  public double getExpirationTime ()          {return T;}
  public double getRiskFreeRate ()            {return r;}
  public double getContinuousDividendYield () {return q;}
  public double getVolatility ()              {return ??;}

  public void print ()
  {
    if(type == true)       { Console.WriteLine("This option is a call"); }
    else if(type == false) { Console.WriteLine("This option is a put"); }
    if(kind == true)       { Console.WriteLine("This option is an American."); }
    else if(kind == false) { Console.WriteLine("This option is an European."); }
    Console.WriteLine("The Strike Price is " + K);
    Console.WriteLine("The Expiration Time is " + T + " years");
    Console.WriteLine("The Risk Free Rate is " + r * 100 + "%");
    Console.WriteLine("The Continuous Dividend Yield is " + q * 100 + "%");
    Console.WriteLine("The Implied Volatility is " + ?? * 100 + "%");
  }

  public double Valuation(double S0, double t0, int n)
  {
      double dt = (T - t0) / (double)n;
      double u = Math.Exp(?? * Math.Sqrt(dt));
      double d = Math.Exp(- ?? * Math.Sqrt(dt));
      double df = Math.Exp(- r * dt);
      double p_rn = (Math.Exp(dt * (r - q)) - d) / (u - d);
      double q_rn = (u - Math.Exp(dt * (r - q))) / (u - d);
      
      double[][] S = null;
      double[][] V = null;
      S = new double[n + 1][];
      V = new double[n + 1][];
      for (int i = 0; i <= n; i++)
      {
          S[i] = new double[i + 1];
          V[i] = new double[i + 1];
      }
      
      S[0][0] = S0;
      for (int i = 1; i <= n; i++)
      {
          for (int j = 0; j <= i; j++)
          {
              S[i][j] = S0 * Math.Pow(u, (i - j)) * Math.Pow(d, j);
          } 
      }
      
      for (int j = 0; j <= n; j++)
      {
          if(type == true)       { V[n][j] = Math.Max((S[n][j] - K), 0); }
          else if(type == false) { V[n][j] = Math.Max((K - S[n][j]), 0); }
      }
      
      for (int i = n - 1; i >= 0; i--)
      {
          for (int j = 0; j <= i; j++)
          {   // discounted expectation
              double Vtmp = df * (p_rn * V[i + 1][j] + q_rn * V[i + 1][j + 1]);
              //Overwrite if it's American Option.
              if(kind == true)
              {
                  if(type == true)  { Vtmp = Math.Max(Vtmp, Math.Max((S[i][j] - K), 0)); }
                  if(type == false) { Vtmp = Math.Max(Vtmp, Math.Max((K - S[i][j]), 0)); }
              }
              
              V[i][j] = Vtmp;
          }
      }
      
      return V[0][0];
  }
  
  public double Del(double S0, double t0, int n)
  {
      double delta = (Valuation((S0 + S0 / 1000), t0, n) - Valuation((S0 - S0 / 1000), t0, n)) / (2 * S0 / 1000);
      
      return delta;
  }
  
  public double Gam(double S0, double t0, int n)
  {
      double gamma = (Valuation((S0 + S0 / 1000), t0, n) + Valuation((S0 - S0 / 1000), t0, n) - 2 * Valuation(S0, t0, n)) / Math.Pow(S0 / 1000, 2);
      
      return gamma;
  }
  
  public double Veg(double S0, double t0, int n)
  {
      double sigma1 = ?? + ?? / 1000;
      double sigma2 = ?? - ?? / 1000;
      
      double ??tmp = ??;
      
      ?? = sigma1;
      double vvega1 = Valuation(S0, t0, n);
      
      ?? = sigma2;
      double vvega2 = Valuation(S0, t0, n);
      
      double vega = (vvega1 - vvega2) / (2 * ??tmp / 1000);
      
      ?? = ??tmp;
      
      return vega;
  }
  
  public double Rhoo(double S0, double t0, int n)
  {
      double r1 = r + r / 1000;
      double r2 = r - r / 1000;
      
      double rtmp = r;
      
      r = r1;
      double vrho1 = Valuation(S0, t0, n);
      
      r = r2;
      double vrho2 = Valuation(S0, t0, n);
      
      double rho = (vrho1 - vrho2) / (2 * rtmp / 1000);
      
      r = rtmp;
      
      return rho;
  }
  
  public double Thet(double S0, double t0, int n)
  {
      double theta = (Valuation(S0, (t0 + (T - t0) / (100 * (double)n)), n) - Valuation(S0, t0, n)) / ((T - t0) / (100 * (double)n));
      
      return theta;
  }
  
  public double impliedVolatility(double S0, double t0, double TargetPrice, double tolerance, int n)
  {
      for(int i = 0; ;i++)
      {
          double ivalue = Valuation(S0, t0, n);
          double ivega = Veg(S0, t0, n);
          
          if (Math.Abs(ivalue - TargetPrice) < tolerance) { break; }
          
          ?? = ?? - (ivalue - TargetPrice) / ivega;
      }
      return (?? * 100);
  }
}
