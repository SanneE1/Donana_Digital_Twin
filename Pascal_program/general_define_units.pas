unit general_define_units;

{$mode ObjFPC}{$H+}

interface

uses
  Classes, SysUtils;

type           {here you declare the data structure for you individuals}
  Array2DInteger = array of array of integer;
  Array3DInteger = array of array of array of integer;

  MapOfLists = array of array of Tlist;

  PRabbit = ^RabbitAgent;

  RabbitAgent = record
    sex: string[1];
    Age: integer;      // months
    Pregnant: boolean;  // Individual pregnant?
    Lactating: boolean;
  end;


 TRHDCellData = record
    NextOutbreakMonth: Integer;      // When next outbreak occurs
    OutbreakFrequency: Integer;      // 24-36 months between outbreaks
    OutbreakMortality: Real;         // 0.15-0.25 mortality rate
  end;



var
  {Simulation info}
  current_year, n_years, start_year, end_year, month, day: integer;
  RabbitMap: Array2DInteger;
  max_pop_size: integer;
  sum_pop_size: array[1..100] of integer;
  each_pop_sizes: array of array of integer;
  output_dir, habitat_folder: string;
  create_maps_3month_maps: boolean = False;
  create_2yr_month_maps: boolean = False;
  rabbit_census_maps: boolean = True;
  create_all_month_maps: boolean = False;

  {General variables}
  Mapdimx, Mapdimy: integer;
  filename: Text;
  xp, yp: integer;

  {Rabbit - population info}
  RabbitPopulationSpatial: MapOfLists;
  RabbitPopulation: TList;
  Rabbit: PRabbit;
  RabbitPopulationSize, RabbitPopulationSpatialSize: longint;
  HabitatMapRabbit: Array2Dinteger;
  paramname_rabbit: string;

  {Rabbit - Climate effects}
  BreedingMonths_xy: Array3DInteger;
  DryMonths_xy: Array3DInteger;
  mapname_rabbit, breedname, dryname, climate_folder: string;
  nbadmonths: byte;

  {Rabbit - Disease mechanism}
  RHDFreq_L, RHDFreq_U: integer;
  RHDMort_L, RHDMort_U: real;
  RHDData: array of array of TRHDCellData;
  tot_m: longint;

  {Rabbit - Demography}
  Kdensity, KCapacity, kCapacity_low, kCapacity_high: real;
  R_lambda, R_dens_opt, R_sigma: real;                            // Dispersal parameters
  R_max_age, R_young_age, R_adult_age, R_old_age: integer;      // Ages in months
  R_disp_age: integer;
  R_MeanLitSize, R_SdLitSize: real;
  R_r_int, R_r_dens_effect, R_r_second_effect, R_r_later_effect: real;
  R_MortP_at_month_old, R_s_int, R_s_extra_juv, R_s_food_effect, R_s_dens_effect: real;

  {Missceleneous // or however you spell that}
   a, i, b, taskID:integer;
   output_maps: string;


const
    days_in_month: array[1..12] of integer = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31);

    L_max_steps = 100;

  disp_x: array[0..48] of integer = (-3, -2, -1, 0, 1, 2, 3,
                                       -3, -2, -1, 0, 1, 2, 3,
                                       -3, -2, -1, 0, 1, 2, 3,
                                       -3, -2, -1, 0, 1, 2, 3,
                                       -3, -2, -1, 0, 1, 2, 3,
                                       -3, -2, -1, 0, 1, 2, 3,
                                       -3, -2, -1, 0, 1, 2, 3);
    disp_y: array[0..48] of integer = (3, 3, 3, 3, 3, 3, 3,
                                       2, 2, 2, 2, 2, 2, 2,
                                       1, 1, 1, 1, 1, 1, 1,
                                       0, 0, 0, 0, 0, 0, 0,
                                      -1, -1, -1, -1, -1, -1, -1,
                                      -2, -2, -2, -2, -2, -2, -2,
                                      -3, -3, -3, -3, -3, -3, -3);
    disp_dist: array[0..48] of integer = (3, 3, 3, 3, 3, 3, 3,
                                          3, 2, 2, 2, 2, 2, 3,
                                          3, 2, 1, 1, 1, 2, 3,
                                          3, 2, 1, 0, 1, 2, 3,
                                          3, 2, 1, 1, 1, 2, 3,
                                          3, 2, 2, 2, 2, 2, 3,
                                          3, 3, 3, 3, 3, 3, 3);

implementation


end.

