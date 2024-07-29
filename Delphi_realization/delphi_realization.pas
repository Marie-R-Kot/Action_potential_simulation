unit Unit1;

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, StdCtrls, DDSpinEdit, Spin, ExtCtrls, TeeProcs, TeEngine, Chart,
  Series, Math, Buttons;

type
  TForm1 = class(TForm)
    Chart1: TChart;
    Series1: TLineSeries;
    Button1: TButton;
    Chart2: TChart;
    Series2: TLineSeries;
    Series3: TLineSeries;
    Chart3: TChart;
    Series4: TLineSeries;
    Series5: TLineSeries;
    DDSpinEdit1: TDDSpinEdit;
    Label1: TLabel;
    Chart4: TChart;
    Series6: TLineSeries;
    Chart5: TChart;
    Series7: TLineSeries;
    Chart6: TChart;
    Series8: TLineSeries;
    Label2: TLabel;
    Label3: TLabel;
    procedure Button1Click(Sender: TObject);
  private
    procedure PerformGaussElimination(p: Integer);
    procedure ExecuteCableSimulation;
  public
  end;

var
  Form1: TForm1;
  a, b, d, f: array[1..9991] of Double;

implementation

{$R *.dfm}

procedure TForm1.PerformGaussElimination(p: Integer);
var
  i: Integer;
begin
  d[1] := d[1] / b[1];
  f[1] := f[1] / b[1];
  b[1] := 1.0;

  for i := 2 to p do
  begin
    b[i] := b[i] - a[i] * d[i - 1];
    f[i] := (f[i] - a[i] * f[i - 1]) / b[i];
    d[i] := d[i] / b[i];
    b[i] := 1.0;
  end;

  for i := 1 to p - 1 do
  begin
    f[p - i] := f[p - i] - f[p - i + 1] * d[p - i];
  end;
end;

procedure TForm1.ExecuteCableSimulation;
var
  C, t, gl, Vl, Ina, Ik, gna, gk, Vna, Vk, Vus, gnap, gkp, dt: Double;
  a1, a2, a3, b1, b3, c1, c3, d1, d2, d3, e1, e2, e3, f1, f2, f3: Double;
  alpham, betam, alphah, betah, alphan, betan, Iel, nu, ri, dx, Rii, Rm, L: Double;
  Length, TimeQ, alpham1, betam1, alphah1, betah1, m1, h1, Ica, gcap, ms1, hs1: Double;
  Vca, Ica1, Integral, m2, alpham2, betam2, ms2, V12, V12na, z, z1: Double;
  time: array[1..300000] of Double;
  R, x, V, Vs, m, h, n, ms, hs, ns, s, u, Vstat: array[1..9991] of Double;
  i, j, p, q, point: Integer;
  outfile: TextFile;

const
  dxValue = 0.01;
  dtValue = 0.01;
  riValue = 0.01;
  LValue = 0.01;
  InitialIntegral = 0.0;
  StartTime = 70;
  EndTime = 75;

begin
  // Initialize variables
  dx := dxValue;
  dt := dtValue;
  ri := riValue;
  L := LValue;
  // Set initial conditions and constants
    Vna:=50;      {mV}
    Vk:=-77;       {mV}
    gnap:=120;     {mS/cm2}
    gkp:=36;       {mS/cm2}
    C:=1;          {mkF/cm2}
    gl:=0.3;       {mS/cm2}
    Vl:=-50;//-50;        {mV}
    Vus:={-65;//}-40;
    Vca:=135;
    V12:=-55;
    V12na:=-50;
    gcap:=4;       {mS/cm2}
    a1:=0.1;
    a2:=25;
    a3:=10;
    b1:=4;
    b3:=18;
    c1:=0.07;
    c3:=20;
    d1:=1;
    d2:=30;
    d3:=10;
    e1:=0.01;
    e2:=10;
    e3:=10;
    f1:=0.125;
    f2:=0;
    f3:=80;
    dt:=0.005;       {ms}
    Iel:=0;         {mkA/cm2}
    Rii:=0.1;      {kOm*cm}
    Rm:=10;        {kOm*cm2}
    Length:=0.1;   {cm}
    L:=2e-3;
    TimeQ:=100;
    q:=3160;
    z:=5;
    z1:=3.4;
    Form1.Series2.Clear;
    Form1.Series3.Clear;
    Form1.Series6.Clear;
    Form1.Series7.Clear;
    Form1.Series8.Clear;
    dx:=Length/(q-1);       {cm}

    point:=trunc(Length/dx)-trunc(L/dx);

    for i:=1 to q  do begin
      x[i]:=(i-1)*Length/(q-1);    {cm}
      R[i]:=1e-4;        {cm}
      V[i]:=0;              {mV}
      Vs[i]:=-61.5;  // }-55;      {mV}
      Vstat[i]:=Vs[i];
	  s[i]:=0;              {mkA/cm2}
      if (i>point-trunc(L/dx)) and (i<point+trunc(L/dx)) then s[i]:=Form1.DDSpinEdit1.Value;//7.646395;  {mS/cm2}
	  if (i>point+trunc(L/dx)) then Vs[i]:= -891373*x[i]*x[i]*x[i] + 304171*x[i]*x[i] - 34101*x[i] + 1210.4;//-1E+09x4 + 5E+08x3 - 8E+07x2 + 5E+06x - 119455-1E+09*x[i]*x[i]*x[i]*x[i] + 6E+08*x[i]*x[i]*x[i] - 8E+07*x[i]*x[i] + 5E+06*x[i]  - 123013;
      if (i<point-trunc(L/dx)) then Vs[i]:=981824*x[i]*x[i]*x[i]*x[i] - 118408*x[i]*x[i]*x[i] + 5398.2*x[i]*x[i] - 75.642*x[i] - 60.999;
      if (i>point-trunc(L/dx)) and (i<point+trunc(L/dx)) then Vs[i]:= -228061*x[i]*x[i] + 41184*x[i] - 1905.1 ;

      u[i]:=0;              {mkA/cm2}
      a[i]:=0;
      b[i]:=0;
      d[i]:=0;
      f[i]:=0;

      alpham:=0.1*(Vs[i]+35)/(1.0-exp(-0.1*(Vs[i]+35))+1e-9);
      betam:=4.0*exp(-0.0556*(Vs[i]- V12na));
      alphah:=0.07*exp(-(Vs[i]- V12na)/20);
      betah:=1.0/(1+exp(-0.1*(Vs[i]- V12na-30))+1e-9);
      alphan:=0.01*(Vs[i]+50)/(1-exp(-0.1*(Vs[i]+50))+1e-9);
      betan:=0.125*exp(-(Vs[i]+60)/80);

      ms[i]:=alpham/(alpham+betam);
      hs[i]:=alphah/(alphah+betah);
      ns[i]:=alphan/(alphan+betan);

    end;

    ri:=Rii/R[1]/R[1];     {kOm/cm}
    for j:=1 to trunc(TimeQ/dt) do begin
        time[j]:=dt*(j-1);
    end;

  // Main simulation loop
  for j := 1 to TimeQ do
  begin
	for i:=1 to q do begin
		point:=q;
		alpham:=0.1*(Vs[i]+35)/(1.0-exp(-0.1*(Vs[i]+35))+1e-9);
		betam:=4.0*exp(-0.0556*(Vs[i]- V12na));
		alphah:=0.07*exp(-(Vs[i]- V12na)/20);
		betah:=1.0/(1+exp(-0.1*(Vs[i]- V12na-30))+1e-9);

		alphan:=0.01*(Vs[i]+50)/(1-exp(-0.1*(Vs[i]+50))+1e-9);
		betan:=0.125*exp(-(Vs[i]+60)/80);

		m[i]:=(alpham+ms[i]/dt)/(alpham+betam+1/dt);
		n[i]:=(alphan+ns[i]/dt)/(alphan+betan+1/dt);
		h[i]:=(alphah+hs[i]/dt)/(alphah+betah+1/dt);

		gna:=gnap*m[i]*m[i]*m[i]*h[i];
		gk:=gkp*n[i]*n[i]*n[i]*n[i];

		a[i]:=1.0/R[i]/2/ri/dx/dx;
		b[i]:=-1.0/R[i]/ri/dx/dx-C/dt-gna-gk-gl-s[i];
		d[i]:=1.0/R[i]/2/ri/dx/dx;
		f[i]:=-C*Vs[i]/dt-gna*Vna-gk*Vk-gl*Vl-s[i]*Vus-u[i];

	end;

		a[1]:=0;  b[1]:=1/ri/dx;  d[1]:=-1/ri/dx;  f[1]:=0; {mkA}
		a[q]:=-1/ri/dx;  b[q]:=1/ri/dx;  d[q]:=0;  f[q]:=0;
		if (dt*j>70) and (dt*j<75) then f[1]:=0.0002;
	  
    // Perform GAUSS elimination for this step
    PerformGaussElimination(q);

    // Update voltage and other parameters
    for i:=1 to q do begin
        V[i]:=f[i];
    end;

    if ( time[j]<=30 ) then begin
		for i:=1 to q do begin
          Vstat[i]:=V[i];
        end;
    end;

    // Output results at specific intervals
    if (j mod 100 = 0) then
    begin
      Form1.Series1.Clear;
      Form1.Series4.Clear;
      Form1.Series5.Clear;
      for i := 1 to q do
      begin
        Form1.Series1.AddXY(dx * (i - 1), V[i]);
        Form1.Series4.AddXY(dx * (i - 1), h[i]);
      end;
      Form1.Series2.AddXY(dt * j, V[Trunc(q / 5)]);
      Form1.Series3.AddXY(dt * j, V[point]);
      Form1.Series6.AddXY(dt * j, Ica);
      Form1.Series7.AddXY(dt * j, m1);
      Form1.Series8.AddXY(dt * j, h1);
      Application.ProcessMessages;
    end;

    // Update state variables for the next iteration
    for i := 1 to q do
    begin
      Vs[i] := V[i];
      ms[i] := m[i];
      ns[i] := n[i];
      hs[i] := h[i];
    end;

    // Update labels with current values
    Form1.Label2.Caption := FloatToStr(V[q]);
    Form1.Label3.Caption := FloatToStr(Integral);
  end;
end;

procedure TForm1.Button1Click(Sender: TObject);
begin
  ExecuteCableSimulation;
end;

end.
