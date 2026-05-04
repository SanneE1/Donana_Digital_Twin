unit rabbit_vital_rates;

{$mode ObjFPC}{$H+}

interface

uses
  Classes, SysUtils, math,
  general_define_units;


procedure Reproduction;
procedure Survival;
procedure Dispersal;

implementation



procedure Reproduction;
var
  a, j, k, minlitsize, maxlitsize: integer;
  rep_p: real;

begin
  with RabbitPopulation do
  begin
    RabbitPopulationSize := RabbitPopulation.Count;
    if (RabbitPopulationSize > 1) then
      for a := 0 to RabbitPopulationSize - 1 do
      begin
        Rabbit := Items[a];

        if Rabbit^.Sex = 'f' then
        {If Rabbit is pregnant and full term; calculate nr of pups and add to population}
        if Rabbit^.Pregnant = True then
        begin
          k:=Round(RandG(R_MeanLitSize,R_SdLitSize));
          minlitsize:=Round(R_MeanLitSize-2*R_SdLitSize);  // Not exactly sure why these 4 lines were in Zulima's code, but I'm keeping it just in case
          maxlitsize:=Round(R_MeanLitSize+2*R_SdLitSize);
          If k < minlitsize then k:= minlitsize;
          If k > maxlitsize then k:= maxlitsize;

          for j:=1 to k do //litter_size_number
            begin
              New(Rabbit);

              Rabbit^.age:=0;
              if random<0.5 then Rabbit^.sex:='f' else Rabbit^.sex:='m';

              Rabbit^.Pregnant:=False;

              RabbitPopulation.Add(Rabbit);
            end;
          Rabbit^.Pregnant:=False
        end
        else
        {If not pregnant, determine if it happens now}
        if (Rabbit^.Age > R_juv_age) then
        begin
        rep_p := 0;

        if (Rabbit^.Age < R_sub_adult_age) then
        begin
           rep_p:= exp(R_r_int - R_r_dens_effect * Kdensity) /
           (1 + exp(R_r_int - R_r_dens_effect * Kdensity)) ;
        end
        else if Rabbit^.Age < R_adult_age then
        begin
          rep_p:= exp(R_r_int + R_r_second_effect - R_r_dens_effect * Kdensity) /
           (1 + exp(R_r_int + R_r_second_effect - R_r_dens_effect * Kdensity));
        end
        else
        begin
          rep_p:= exp(R_r_int + R_r_later_effect - R_r_dens_effect * Kdensity) /
           (1 + exp(R_r_int + R_r_later_effect - R_r_dens_effect * Kdensity));
        end;


        if random < rep_p then
        begin
          Rabbit^.Pregnant:=True;
        end;
      end;
    end;
  end;
end;

procedure Survival;
var
  a: integer;
  surv_nb, surv_juv, surv_adult, surv_p: real;
  die: boolean;

begin
  {Calculate current survival probabilities as a function of density}

        surv_nb:=exp(ln(1-R_MortP_at_month_old));

        surv_juv:=(exp(R_s_int - R_s_extra_juv -R_s_dens_effect*Kdensity)/
                 (1+exp(R_s_int - R_s_extra_juv -R_s_dens_effect*Kdensity)))*
                 (1-exp(-R_s_food_effect/(nbadmonths+0.001)));

        surv_adult:=(exp(R_s_int-R_s_dens_effect*Kdensity)/
                 (1+exp(R_s_int-R_s_dens_effect*Kdensity)))*
                 (1-exp(-R_s_food_effect/(nbadmonths+0.001)));

  {Determine survival of Rabbits}
  with RabbitPopulation do
  begin

    for a := RabbitPopulationSize - 1 downto 0 do  //By going down instead of up, I don't have to continuously correct the index when an Rabbit dies
    begin
      Rabbit := items[a];

      {Assign current survival probabilities}
      if Rabbit^.age = 0 then surv_p := surv_nb
      else if Rabbit^.age <= R_sub_adult_age then surv_p := surv_juv
      else surv_p := surv_adult;

      {Determine fate of Rabbits}
      die := False;
      if Rabbit^.age >= R_max_age then die := True
      else
      if random > surv_p then die := True;
      if die then
        Delete(a);
    end;

  end;
end;

procedure Dispersal;
var
iy, ix, b, c, dens: integer;
p_not_stand, p_stand, p_cumsum: array of real;
habitable: array of boolean;
moved: boolean;
p_total, p_sum, p_disp: real;
begin

  // Loop over the spatial population and find out how many Rabbits are dispersing
  for iy := 1 to MapdimY do
  begin
    for ix := 1 to MapdimX do
    begin

      if RabbitPopulationSpatial[ix, iy] = nil then Continue;

      RabbitPopulation := RabbitPopulationSpatial[ix, iy];
      RabbitPopulationSize := RabbitPopulation.Count;

      setLength(p_not_stand, 49);
      setLength(p_stand, 49);
      setLength(p_cumsum, 49);

      //setLength(habitable, 49);
      //FillChar(habitable[0], Length(habitable) * SizeOf(habitable), True);

      p_total := 0;
      p_sum := 0;

     for c:=0 to Length(disp_x) - 1 do
     begin
      if (ix + disp_x[c] < 1) or (iy + disp_y[c] < 1) or
      (ix + disp_x[c] > MapDimX) or (iy + disp_y[c] > MapDimY) or
      (HabitatMapRabbit[ix + disp_x[c], iy + disp_y[c]] = 0) then
      begin
      p_not_stand[c] := 0;
      //habitable[c] := false;
      end
      else
      begin
       if RabbitPopulationSpatial[ix + disp_x[c], iy + disp_y[c]] = nil then
       dens := 0
       else
       dens := RabbitPopulationSpatial[ix + disp_x[c], iy + disp_y[c]].Count;
       //if dens = nil then dens := 0;
       p_not_stand[c] := exp(-R_lambda * ((dens) - R_dens_opt)**2) / ((1 + disp_dist[c])**R_sigma);
      end;
     end;

     for c:=0 to Length(disp_x) - 1 do
      p_total:= p_total + p_not_stand[c];

     for c:=0 to Length(disp_x) - 1 do
     begin
     p_stand[c] := p_not_stand[c] / p_total;
     p_sum:= p_sum + p_stand[c];
     p_cumsum[c] := p_sum;
     end;

     with RabbitPopulation do
          for b := RabbitPopulationSize - 1 downto 0 do
        begin
          Rabbit := Items[b];

          if Rabbit^.Age = R_disp_age then
          begin
          moved := false;

          p_disp := random;

          for c:=0 to Length(disp_x) - 1 do
          if {(habitable[c]) and }(p_disp < p_cumsum[c]) then
          begin
           moved := true;
           if (ix + disp_x[c] < 1) or (iy + disp_y[c] < 1) or (ix + disp_x[c] > MapDimX) or (iy + disp_y[c] > MapDimY) then
          Exit;

          {Remove Rabbit from current population}
          Delete(b);

          {Moving rabbit to new location}
          if RabbitPopulationSpatial[ix + disp_x[c], iy + disp_y[c]] = nil then
          RabbitPopulationSpatial[ix + disp_x[c], iy + disp_y[c]] := TList.Create;

          RabbitPopulationSpatial[ix + disp_x[c], iy + disp_y[c]].Add(Rabbit);

          {if check_disp_count = 0 then
          setLength(check_disp, 20000);
          if check_disp_count < 20000 then
          begin
          check_disp[check_disp_count] := c;
          Inc(check_disp_count);
          end;}

          {stop loop}
          Break;

          end;

          if moved = false then
          Exit;

          end;
        end;
    end;
  end;
end;

end.

