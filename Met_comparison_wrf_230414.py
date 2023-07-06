# -*- coding: utf-8 -*-
# Preceding script: E:\JGY\Code\Important_Code_data\Emis_Process\Emis_allocation_ABaCAS_9k.py (
# Modified version from Matlab (E:\JGY\Code\Important_Code_data\WRF_Validation\WRF_Simulation_Validation_230401.m)
#######################################################################
# Data Format of isd ground met observations:
#######################################################################
# Field 1: Pos 1-4, Length 4: Observation Year
# Year of observation, rounded to nearest whole hour
#
# Field 2: Pos 6-7, Length 2: Observation Month
# Month of observation, rounded to nearest whole hour
#
# Field 3: Pos 9-11, Length 2: Observation Day
# Day of observation, rounded to nearest whole hour
#
# Field 4: Pos 12-13, Length 2: Observation Hour
# Hour of observation, rounded to nearest whole hour
#
# Field 5: Pos 14-19, Length 6:  Air Temperature
# The temperature of the air
# UNITS:  Degrees Celsius
# SCALING FACTOR: 10
# MISSING VALUE: -9999
#
# Field 6: Pos 20-24, Length 6: Dew Point Temperature The temperature to which a given parcel of air must be cooled
# at constant pressure and water vapor content in order for saturation to occur. UNITS: Degrees Celsius SCALING
# FACTOR: 10 MISSING VALUE: -9999
#
# Field 7: Pos 26-31, Length 6: Sea Level Pressure
# The air pressure relative to Mean Sea Level (MSL).
# UNITS:  Hectopascals
# SCALING FACTOR: 10
# MISSING VALUE: -9999
#
# Field 8: Pos 32-37, Length 6: Wind Direction
# The angle, measured in a clockwise direction, between true north and the direction from which the wind is blowing.
# UNITS: Angular Degrees
# SCALING FACTOR: 1
# MISSING VALUE: -9999
# *NOTE:  Wind direction for calm winds is coded as 0.
#
# Field 9: Pos 38 - 43, Length 6: Wind Speed Rate
# The rate of horizontal travel of air past a fixed point.
# UNITS: meters per second
# SCALING FACTOR: 10
#######################################################################
"""
本脚本WRF结果验证思路：
1. 对于模拟值提取逻辑：
   loop through 站点列表中的站点：
    if 当前站点位于domain:
        从METCRO2D文件中提取所有观测要素（T2，WS，WD，P0,），检查参数的单位，做好必要的转换；
        从Met_grid文件中提取RH；
        将模拟结果、UTC年月日时（分）保存至数组M；
    else:
        下一个站点

   保存数组M至文件Model_SiteID_Para_Year。

2. 站点值提取：

   # 检查中间文件是否已存在
   if Obs中间文件已存在:
    跳过第2步
   else:
    loop through 站点数据文件夹中的文件:
        读取当前文件txt, 将关心的变量（）以及UTC时间存至数组Obs，检查参数的单位，做好必要的转换；

   保存Obs至文件Obs_SiteID_Year

3. Site-Grid 关联与比较
   采用Python相应模块进行数据比较。
# % Add linear correlation coefficients (R).

更新日志 Update History
Apr 14, 2023
directly reading the wrfout file.

"""

import time
import numpy as np
import matplotlib
from netCDF4 import Dataset
from datetime import datetime, date, timedelta
import pandas as pd
from scipy.io import savemat, loadmat
import math
import glob


def read_ascii(file):
    """
    读取ASCII文件
    :param file: 输入是f.readlines(),一个列表文件
    :return: 还是一个列表文件
    :功能：防止文件中间或者末尾有空行
    本函数参考自 https://blog.csdn.net/qq_44907989/article/details/125907641
    """

    file_tep = []
    for j in range(len(file)):
        if len(file[j].split()) != 0:
            file_tep.append(file[j])
    file = file_tep
    return file


"""
%% Function rh_cal for Relative Humidity (in percentage)
% p: surface pressure (in Pa)
% q: water vapor mixing ratio (kg/kg)
% t: temperature at 2 meter a.g. (Kelvin)
"""


def rh_cal(p, q, t):
    svp1 = 611.2
    svp2 = 17.67
    svp3 = 29.65
    svpt0 = 273.15
    eps = 0.622
    rh = 100 * (p * q / (q * (1 - eps) + eps)) / (svp1 * np.exp(svp2 * (t - svpt0) / (t - svp3)))
    rh[rh > 100] = 100
    return rh


"""
Function rh_from_dt_cal calculating the rh based on the dew point and air temperature
dp: dew point temp (Celsius Deg.)
tempa: air temp (Celsius Deg.)
"""


def rh_from_dt_cal(dp, tempa):
    # beta and lamda -  Revised Magnus coefficients as recommended by Alduchov and Eskridge
    # doi.org/10.1175/1520-0450(1996)035<0601:IMFAOS>2.0.CO;2
    # the formula is directly copied from website: https://www.omnicalculator.com/physics/dew-point
    beta = 17.625
    lamda = 243.04  # Unit is Celsius degree
    rh = 100 * np.exp((beta * dp) / (lamda + dp)) / np.exp((beta * tempa) / (lamda + tempa))  # in %
    return rh


if __name__ == '__main__':
    a = time.time()
    # some basic configuration information
    met_grid_ncs_path = r""
    # 站点观测数据文件夹
    site_obs_path = r"E:\CMAQ_run_data\Met_data\china_isd_lite_2017\china_isd_lite_2017"
    # 站点经纬度信息
    site_info_dir = r"E:\CMAQ_run_data\Met_data\_数据说明\isd-history_CN.csv"
    # metcro_2d_path = r"E:\CMAQ_run_data\230412\python_val_test\METCRO2D_171102_d01.nc"
    # wrfout文件
    metcro_2d_path = r"E:\CMAQ_run_data\230412\wrfout_d03_2017-04-01_09_00_00.nc"
    gridcro_2d_pat = metcro_2d_path
    # 输出结果文件（主要包括分站点的验证结果以及所有站点的验证结果，CSV格式）
    output_var_pat = r"E:\CMAQ_run_data\230415\d01_val"

    # Meteorological result Simulation extraction
    # loop through the sites in "站点列表"
    f = open(site_info_dir)
    file = read_ascii(f.readlines())
    rown = len(file)
    # three columns in site_info denote respectively SiteID, Lon, Lat
    # refer to "E:\CMAQ_run_data\Met_data\_数据说明\isd-history.txt" for details
    site_info = np.zeros((rown - 1, 3))
    current_row_number = 0
    for num, value in enumerate(file):
        value = value.split(',')
        if num > 0:
            site_info[current_row_number, 0] = float(value[0])
            if (value[6] == '') | (value[7] == ''):
                continue
            site_info[current_row_number, 2] = float(value[6])
            site_info[current_row_number, 1] = float(value[7])
            current_row_number = current_row_number + 1

    f.close()
    site_info = site_info[~np.all(site_info == 0, axis=1)]

    # extract the lon-lat grid info from GRIDCRO2D

    grid_lonlat = Dataset(gridcro_2d_pat, 'r')
    grid_Metcro = Dataset(metcro_2d_path, 'r')
    lon_grid = grid_lonlat.variables['XLONG'][:]
    lat_grid = grid_lonlat.variables['XLAT'][:]
    lon_grid = lon_grid[0, :, :]
    lat_grid = lat_grid[0, :, :]
    T2 = grid_Metcro.variables['T2'][:]
    TFLAG = grid_Metcro.variables['Times'][:]  # [:, 0, :]
    WSped10_u = grid_Metcro.variables['U10'][:]
    WSped10_v = grid_Metcro.variables['V10'][:]
    # WSped10 =

    Q2 = grid_Metcro.variables['Q2'][:]
    P0 = grid_Metcro.variables['PSFC'][:]
    
    site_num = len(site_info)
    for ss in range(site_num):
        curr_lon = site_info[ss, 1]
        curr_lat = site_info[ss, 2]
        curr_sID = site_info[ss, 0]
        dist = (lon_grid - curr_lon) ** 2 + (lat_grid - curr_lat) ** 2
        min_dist = dist.min()
        if min_dist > 0.09:
            continue
        # print(min_dist)
        ind_2d = np.unravel_index(np.argmin(dist, axis=None), dist.shape)
        T2_timeseries = T2[:, ind_2d[0], ind_2d[1]]  # simplified as "_ts" below
        # WSped10_ts = WSped10[:, 0, ind_2d[0], ind_2d[1]]
        WSped10_u_ts = WSped10_u[:, ind_2d[0], ind_2d[1]]
        WSped10_v_ts = WSped10_v[:, ind_2d[0], ind_2d[1]]
        WSped10_ts = np.sqrt(WSped10_v_ts ** 2 + WSped10_u_ts ** 2)
        # 注意使用arctan2 而不是arctan (Using np.arctan2 instead of arctan)
        WDire10_mathe = np.arctan2(WSped10_v_ts,
                                   WSped10_u_ts) * 180 / np.pi  # np.arctan(WSped10_v_ts / WSped10_u_ts)  #
        # grid_Metcro.variables['WDIR10'][:]
        WDire10_ts_1 = 270 - WDire10_mathe
        WDire10_ts = np.mod(WDire10_ts_1, 360)
        # WDire10_ts_real = np.mod(WDire10_ts, 360)
        # WDire10_ts = WDire10[:, 0, ind_2d[0], ind_2d[1]]
        Q2_ts = Q2[:, ind_2d[0], ind_2d[1]]
        P0_ts = P0[:, ind_2d[0], ind_2d[1]]
        rh_ts = rh_cal(P0_ts, Q2_ts, T2_timeseries)
        # date_list = str(TFLAG[:, 0:19], 'UTF-8')

        # Creating empty series
        # hr_ser = pd.Series([], dtype=int64)
        hr_ser = np.zeros((TFLAG.shape[0], 1), dtype=np.uint32)
        date_ser = np.zeros((TFLAG.shape[0], 1), dtype=np.uint32)
        for ii in range(TFLAG.shape[1]):
            tmutc_list = str(TFLAG[ii, :], 'UTF-8')
            yr = int(tmutc_list[0:4])
            mon = int(tmutc_list[5:7])
            dy = int(tmutc_list[8:10])
            hr = int(tmutc_list[11:13])
            doy = date(yr, mon, dy).timetuple().tm_yday
            hr_1000 = hr * 10000
            hr_ser[ii] = hr_1000
            date_ser[ii] = doy + yr * 1000
            # print(tmutc_list)

        yr = np.floor(date_ser[0] / 1000)
        sim_res_arr = np.zeros((len(date_ser), 9))
        sim_res_arr[:, 0] = curr_sID
        sim_res_arr[:, 1] = date_ser[:, 0]
        sim_res_arr[:, 2] = hr_ser[:, 0]
        sim_res_arr[:, 3] = T2_timeseries - 273.15  # convert to Celsius degree
        sim_res_arr[:, 4] = WSped10_ts
        sim_res_arr[:, 5] = WDire10_ts
        sim_res_arr[:, 6] = Q2_ts
        sim_res_arr[:, 7] = P0_ts / 100
        sim_res_arr[:, 8] = rh_ts
        output_sim_res_path = output_var_pat + '\\Sim_Site_' + str(int(curr_sID)) + '_year_' + str(int(yr)) + '.pkl'
        df_sim_res_arr = pd.DataFrame(sim_res_arr, columns=['SiteID', 'Date', 'Hour(Utc)', 'Temp2m(Celcius)',
                                                            'Wind Speed at 10m (m s-1)',
                                                            'Wind Direction(Deg)', 'Water vapor mixing ratio',
                                                            'surface pressure (hPa)', 'RH (%)'])
        df_sim_res_arr.to_pickle(output_sim_res_path)
        # mdict = {'sim_res_arr': sim_res_arr}
        # savemat(output_sim_res_path, mdict)

    print("Complete the MCIP extraction part~")

    # Start in situ data reading procession!

    # first, check if the observation files have been processed
    # if so, skip this part and go to the comparison part.
    src_str_pkl = output_var_pat + "\\Obs*.pkl"
    files_obspkl = glob.glob(src_str_pkl)
    pkl_file_count = len(files_obspkl)
    # setting of scale factors of ground based observations 比例系数设置（参考数据说明，或本脚本开头数据介绍部分）
    sf_temp = 10
    sf_dewT = 10
    sf_slpr = 10
    sf_wddir = 1.0
    sf_wdspd = 10

    search_file_pattern = site_obs_path + '\\*-2017'
    Met_Obs_files = glob.glob(search_file_pattern)
    met_file_cnt = len(Met_Obs_files)
    if met_file_cnt == 0:
        print("Check the input ground observations! ")
    elif met_file_cnt == pkl_file_count:
        print("All files have been processed. Skip to Comparison part! ")
    else:
        for ii in range(met_file_cnt):
            current_obs_path = Met_Obs_files[ii]
            current_sID = current_obs_path.split('\\')[-1].split('-')[0]
            f = open(current_obs_path)
            file = read_ascii(f.readlines())
            rown = len(file)
            # THE six fields are respectively DOY, Hour (utc), surface temperature, dew point temperature,
            # sea level pressure,  Wind Direction, and Wind Speed Rate
            grd_obs = np.zeros((rown, 7))
            current_row_number = 0
            for num, value in enumerate(file):
                value = value.split(' ')
                value = [i for i in value if i]
                # if '-9999' in value[4:10]:  # 包含缺失值
                #     continue
                year = int(value[0])
                mon = int(value[1])
                dy = int(value[2])
                doy = date(year, mon, dy).timetuple().tm_yday
                grd_obs[current_row_number, 0] = doy + 1000 * year
                grd_obs[current_row_number, 1] = float(value[3]) * 10000  # hour(utc)
                grd_obs[current_row_number, 2] = float(value[4])  # / sf_temp  # air temp
                grd_obs[current_row_number, 3] = float(value[5])  # / sf_dewT  # dew oiubt temp
                grd_obs[current_row_number, 4] = float(value[6])  # / sf_slpr  # mean sea level pressure
                grd_obs[current_row_number, 5] = float(
                    value[7])  # / sf_wddir  # Wind direction for calm winds is coded as 0.
                grd_obs[current_row_number, 6] = float(value[8])  # / sf_wdspd  # Wind Speed Rate
                current_row_number = current_row_number + 1
            grd_obs = grd_obs[~np.all(grd_obs == 0, axis=1)]
            grd_obs[grd_obs == -9999] = np.nan
            grd_obs[:, 2] = grd_obs[:, 2] / sf_temp
            grd_obs[:, 3] = grd_obs[:, 3] / sf_dewT
            grd_obs[:, 4] = grd_obs[:, 4] / sf_slpr
            grd_obs[:, 5] = grd_obs[:, 5] / sf_wddir
            grd_obs[:, 6] = grd_obs[:, 6] / sf_wdspd
            df_grd_obs = pd.DataFrame(grd_obs, columns=['Date', 'Hour(Utc)', 'Air Temp(Celcius)', 'Dew point temp',
                                                        'Sea level pressure (hPa)', 'wind dir', 'wind speed (m s-1)'])
            f.close()
            output_obs_res_path = output_var_pat + '\\Obs_Site_' + current_sID + '_year_' + str(year) + '.pkl'
            df_grd_obs.to_pickle(output_obs_res_path)
            # mdict = {'df_grd_obs': df_grd_obs}
            # savemat(output_obs_res_path, mdict)
        print("!!Complete Ground Observation!!")

    ## Start ground-satellite comparison!!
    search_mat_path = output_var_pat + r"\Obs*.pkl"
    obs_mats = glob.glob(search_mat_path)
    mat_count = len(obs_mats)
    com_res_array_each_site_t = np.zeros((mat_count, 8))
    com_res_array_each_site_Ws = np.zeros((mat_count, 8))
    com_res_array_each_site_Wd = np.zeros((mat_count, 8))
    com_res_array_each_site_rh = np.zeros((mat_count, 8))
    # com_res_array_each_site_RH = np.zeros((mat_count, 8))

    current_site = 0
    # 结果汇总
    df_sim_t = []
    df_obs_t = []
    df_sim_ws = []
    df_obs_ws = []
    df_sim_wd = []
    df_obs_wd = []
    df_obs_rh = []
    df_sim_rh = []
    start_date_str = []
    end_date_str = []
    for ff in range(mat_count):
        current_mat_Obs = obs_mats[ff]
        date_str = current_mat_Obs.split('\\')[-1].split('_')[2]
        src_sim_path = output_var_pat + r"\Sim*" + date_str + '*.pkl'
        sim_files = glob.glob(src_sim_path)
        if len(sim_files) == 0:
            print("Corresponding site file not found！")
            continue

        # sim_data = loadmat(sim_files[0])
        # obs_data = loadmat(current_mat_Obs)
        sim_data = pd.read_pickle(sim_files[0])
        obs_data = pd.read_pickle(current_mat_Obs)
        comparison_met = pd.merge(obs_data, sim_data, on=['Hour(Utc)', 'Date'])

        if len(comparison_met) == 0:
            continue
        ###########################################
        # TEMPERATURE
        ###########################################
        site_obs_t = comparison_met['Air Temp(Celcius)']
        sim_t = comparison_met['Temp2m(Celcius)']

        start_date_str = str(int(comparison_met['Date'][0])) + '_' + str(int(comparison_met['Hour(Utc)'][0]))
        end_date_str = str(int(comparison_met['Date'][len(comparison_met) - 1])) + '_' + str(
            int(comparison_met['Hour(Utc)'][len(comparison_met) - 1]))
        # For some sites observations are missing, skip to the next site.

        if current_site == 0:
            df_sim_t = sim_t
            df_obs_t = site_obs_t
        else:
            df_sim_t = df_sim_t.append(sim_t)
            df_obs_t = df_obs_t.append(site_obs_t)

        r_t = sim_t.corr(site_obs_t)
        rmse_t = ((site_obs_t - sim_t) ** 2).mean() ** .5
        ge_t = (site_obs_t - sim_t).abs().mean()  # mean(abs(jn_sim_temp - jn_obs_temp), 'omitnan')
        bias_t = (site_obs_t - sim_t).mean()
        avg_sim_t = sim_t.mean()
        avg_obs_t = site_obs_t.mean()
        dev_sim_abs = abs(sim_t - avg_sim_t)
        dev_obs_abs = abs(site_obs_t - avg_obs_t)
        ioa_t = 1 - np.nansum((sim_t - site_obs_t) ** 2) / np.nansum((dev_obs_abs + dev_sim_abs) ** 2)
        siteID = sim_data['SiteID'][0]
        com_res_array_each_site_t[current_site, 0:8] = [siteID, avg_sim_t, avg_obs_t, r_t, rmse_t, ge_t, bias_t, ioa_t]
        ###########################################
        # Wind Speed
        ###########################################
        site_obs_ws = comparison_met['wind speed (m s-1)']
        sim_ws = comparison_met['Wind Speed at 10m (m s-1)']
        if current_site == 0:
            df_sim_ws = sim_ws
            df_obs_ws = site_obs_ws
        else:
            df_sim_ws = df_sim_ws.append(sim_ws)
            df_obs_ws = df_obs_ws.append(site_obs_ws)

        r_ws = sim_ws.corr(site_obs_ws)
        rmse_ws = ((site_obs_ws - sim_ws) ** 2).mean() ** .5
        ge_ws = (site_obs_ws - sim_ws).abs().mean()  # mean(abs(jn_sim_temp - jn_obs_temp), 'omitnan')
        bias_ws = (site_obs_ws - sim_ws).mean()
        avg_sim_ws = sim_ws.mean()
        avg_obs_ws = site_obs_ws.mean()
        dev_sim_abs = abs(sim_ws - avg_sim_ws)
        dev_obs_abs = abs(site_obs_ws - avg_obs_ws)
        ioa_ws = 1 - np.nansum((sim_ws - site_obs_ws) ** 2) / np.nansum((dev_obs_abs + dev_sim_abs) ** 2)
        siteID = sim_data['SiteID'][0]
        com_res_array_each_site_Ws[current_site, 0:8] = [siteID, avg_sim_ws, avg_obs_ws, r_ws, rmse_ws, ge_ws, bias_ws,
                                                         ioa_ws]
        ###########################################
        # Wind Direction
        ###########################################
        site_obs_Wd = comparison_met['wind dir']
        sim_Wd = comparison_met['Wind Direction(Deg)']
        if current_site == 0:
            df_sim_wd = sim_Wd
            df_obs_wd = site_obs_Wd
        else:
            df_sim_wd = df_sim_wd.append(sim_Wd)
            df_obs_wd = df_obs_wd.append(site_obs_Wd)

        r_wd = sim_Wd.corr(site_obs_Wd)
        rmse_wd = ((site_obs_Wd - sim_Wd) ** 2).mean() ** .5
        ge_wd = (site_obs_Wd - sim_Wd).abs().mean()  # mean(abs(jn_sim_temp - jn_obs_temp), 'omitnan')
        bias_wd = (site_obs_Wd - sim_Wd).mean()
        avg_sim_wd = sim_Wd.mean()
        avg_obs_wd = site_obs_Wd.mean()
        dev_sim_abs = abs(sim_Wd - avg_sim_wd)
        dev_obs_abs = abs(site_obs_Wd - avg_obs_wd)
        ioa_wd = 1 - np.nansum((sim_Wd - site_obs_Wd) ** 2) / np.nansum((dev_obs_abs + dev_sim_abs) ** 2)
        siteID = sim_data['SiteID'][0]
        com_res_array_each_site_Wd[current_site, 0:8] = [siteID, avg_sim_wd, avg_obs_wd, r_wd, rmse_wd, ge_wd, bias_wd,
                                                         ioa_wd]
        ###########################################
        # Surface Pressure (两者的定义的一致性存疑)
        ###########################################

        # print(sim_data)
        # print(obs_data)
        # Temp 2m
        ###########################################
        # RH
        ###########################################

        # calculating the ground - RH from air temperature and dew point temperature
        dp_obs = comparison_met['Dew point temp']
        site_obs_t = comparison_met['Air Temp(Celcius)']
        sim_rh = comparison_met['RH (%)']
        # both input temperatures must be in Celsius deg. unit!
        rh_obs = rh_from_dt_cal(dp_obs, site_obs_t)
        if current_site == 0:
            df_sim_rh = sim_rh
            df_obs_rh = rh_obs
        else:
            df_sim_rh = df_sim_rh.append(sim_rh)
            df_obs_rh = df_obs_rh.append(rh_obs)

        r_rh = sim_rh.corr(rh_obs)
        rmse_rh = ((rh_obs - sim_rh) ** 2).mean() ** .5
        ge_rh = (rh_obs - sim_rh).abs().mean()  # mean(abs(jn_sim_temp - jn_obs_temp), 'omitnan')
        bias_rh = (rh_obs - sim_rh).mean()
        avg_sim_rh = sim_rh.mean()
        avg_obs_rh = rh_obs.mean()
        dev_sim_abs = abs(sim_rh - avg_sim_rh)
        dev_obs_abs = abs(rh_obs - avg_obs_rh)
        ioa_rh = 1 - np.nansum((sim_rh - rh_obs) ** 2) / np.nansum((dev_obs_abs + dev_sim_abs) ** 2)
        siteID = sim_data['SiteID'][0]
        com_res_array_each_site_rh[current_site, 0:8] = [siteID, avg_sim_rh, avg_obs_rh, r_rh, rmse_rh, ge_rh, bias_rh,
                                                         ioa_rh]

        current_site = current_site + 1
    com_res_array_each_site_t = com_res_array_each_site_t[~np.all(com_res_array_each_site_t == 0, axis=1)]  # 去除全零行
    df_com_res_array_each_site_t = pd.DataFrame(com_res_array_each_site_t,
                                                columns=["Site ID", "SimAvg", "ObsAvg", "LinearCorr", "RMSE_T2",
                                                         "GrossError", "MeanBias", "IOA"])
    out_csv_t_name = output_var_pat + '\\metric_temp2m_site_based_' + start_date_str + '_' + end_date_str + '.csv'
    df_com_res_array_each_site_t.to_csv(out_csv_t_name, index=False)

    com_res_array_each_site_Ws = com_res_array_each_site_Ws[~np.all(com_res_array_each_site_Ws == 0, axis=1)]  # 去除全零行
    df_com_res_array_each_site_Ws = pd.DataFrame(com_res_array_each_site_Ws,
                                                 columns=["Site ID", "SimAvg", "ObsAvg", "LinearCorr", "RMSE_Ws10",
                                                          "GrossError", "MeanBias", "IOA"])
    out_csv_Ws_name = output_var_pat + '\\metric_Windsp10m_site_based' + start_date_str + '_' + end_date_str + '.csv'
    df_com_res_array_each_site_Ws.to_csv(out_csv_Ws_name, index=False)

    com_res_array_each_site_Wd = com_res_array_each_site_Wd[~np.all(com_res_array_each_site_Wd == 0, axis=1)]  # 去除全零行
    df_com_res_array_each_site_Wd = pd.DataFrame(com_res_array_each_site_Wd,
                                                 columns=["Site ID", "SimAvg", "ObsAvg", "LinearCorr", "RMSE_Wd",
                                                          "GrossError", "MeanBias", "IOA"])
    out_csv_Wd_name = output_var_pat + '\\metric_WdDir10_site_based' + start_date_str + '_' + end_date_str + '.csv'
    df_com_res_array_each_site_Wd.to_csv(out_csv_Wd_name, index=False)

    com_res_array_each_site_rh = com_res_array_each_site_rh[~np.all(com_res_array_each_site_rh == 0, axis=1)]  # 去除全零行
    df_com_res_array_each_site_rh = pd.DataFrame(com_res_array_each_site_rh,
                                                 columns=["Site ID", "SimAvg", "ObsAvg", "LinearCorr", "RMSE_rh",
                                                          "GrossError", "MeanBias", "IOA"])
    out_csv_rh_name = output_var_pat + '\\metric_RH_site_based' + start_date_str + '_' + end_date_str + '.csv'
    df_com_res_array_each_site_rh.to_csv(out_csv_rh_name, index=False)
    # overall evaluation
    com_all_sites_each_indicator = np.zeros((7, 4))  # each column denotes different fields
    # Temp2
    r_t = df_sim_t.corr(df_obs_t)
    rmse_t = ((df_obs_t - df_sim_t) ** 2).mean() ** .5
    ge_t = (df_obs_t - df_sim_t).abs().mean()  # mean(abs(jn_sim_temp - jn_obs_temp), 'omitnan')
    bias_t = (df_obs_t - df_sim_t).mean()
    avg_sim_t = df_sim_t.mean()
    avg_obs_t = df_obs_t.mean()
    dev_sim_abs = abs(df_sim_t - avg_sim_t)
    dev_obs_abs = abs(df_obs_t - avg_obs_t)
    ioa_t = 1 - np.nansum((df_sim_t - df_obs_t) ** 2) / np.nansum((dev_obs_abs + dev_sim_abs) ** 2)
    com_all_sites_each_indicator[:, 0] = [avg_sim_t, avg_obs_t, r_t, rmse_t, ge_t, bias_t, ioa_t]

    # Wspeed10
    site_obs_ws = df_obs_ws
    sim_ws = df_sim_ws
    r_ws = sim_ws.corr(site_obs_ws)
    rmse_ws = ((site_obs_ws - sim_ws) ** 2).mean() ** .5
    ge_ws = (site_obs_ws - sim_ws).abs().mean()  # mean(abs(jn_sim_temp - jn_obs_temp), 'omitnan')
    bias_ws = (site_obs_ws - sim_ws).mean()
    avg_sim_ws = sim_ws.mean()
    avg_obs_ws = site_obs_ws.mean()
    dev_sim_abs = abs(sim_ws - avg_sim_ws)
    dev_obs_abs = abs(site_obs_ws - avg_obs_ws)
    ioa_ws = 1 - np.nansum((sim_ws - site_obs_ws) ** 2) / np.nansum((dev_obs_abs + dev_sim_abs) ** 2)
    com_all_sites_each_indicator[:, 1] = [avg_sim_ws, avg_obs_ws, r_ws, rmse_ws, ge_ws, bias_ws, ioa_ws]

    # Wdir10
    site_obs_Wd = df_obs_wd
    sim_Wd = df_sim_wd
    r_wd = sim_Wd.corr(site_obs_Wd)
    rmse_wd = ((site_obs_Wd - sim_Wd) ** 2).mean() ** .5
    ge_wd = (site_obs_Wd - sim_Wd).abs().mean()  # mean(abs(jn_sim_temp - jn_obs_temp), 'omitnan')
    bias_wd = (site_obs_Wd - sim_Wd).mean()
    avg_sim_wd = sim_Wd.mean()
    avg_obs_wd = site_obs_Wd.mean()
    dev_sim_abs = abs(sim_Wd - avg_sim_wd)
    dev_obs_abs = abs(site_obs_Wd - avg_obs_wd)
    ioa_wd = 1 - np.nansum((sim_Wd - site_obs_Wd) ** 2) / np.nansum((dev_obs_abs + dev_sim_abs) ** 2)
    com_all_sites_each_indicator[:, 2] = [avg_sim_wd, avg_obs_wd, r_wd, rmse_wd, ge_wd, bias_wd, ioa_wd]

    # RH
    sim_rh = df_sim_rh
    rh_obs = df_obs_rh

    r_rh = sim_rh.corr(rh_obs)
    rmse_rh = ((rh_obs - sim_rh) ** 2).mean() ** .5
    ge_rh = (rh_obs - sim_rh).abs().mean()  # mean(abs(jn_sim_temp - jn_obs_temp), 'omitnan')
    bias_rh = (rh_obs - sim_rh).mean()
    avg_sim_rh = sim_rh.mean()
    avg_obs_rh = rh_obs.mean()
    dev_sim_abs = abs(sim_rh - avg_sim_rh)
    dev_obs_abs = abs(rh_obs - avg_obs_rh)
    ioa_rh = 1 - np.nansum((sim_rh - rh_obs) ** 2) / np.nansum((dev_obs_abs + dev_sim_abs) ** 2)
    com_all_sites_each_indicator[:, 3] = [avg_sim_rh, avg_obs_rh, r_rh, rmse_rh, ge_rh, bias_rh, ioa_rh]

    df_com_all_sites_each_indicator = pd.DataFrame(com_all_sites_each_indicator,
                                                   columns=["Temp2", "WindSpeed@10m", "WindDir@10m", "RH(%)"],
                                                   index=["avg_sim", "avg_obs", "R", "rmse", "GE", "MB", "IOA"])
    out_csv_total_name = output_var_pat + '\\Allsites_all_indic_d01_' + start_date_str + '_' + end_date_str + '.csv'
    df_com_all_sites_each_indicator.to_csv(out_csv_total_name)

    b = time.time()
    print("总时间:", b - a)
    # day_of_year = date(2017, 1, 31).timetuple().tm_yday
    # print(day_of_year)
