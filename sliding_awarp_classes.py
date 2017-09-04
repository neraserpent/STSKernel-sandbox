class SlidingAwarp():

    def __init__(self):
        pass
    
    def _AWarp(self, x, y):
        import numpy as np
        n = len(x)
        m = len(y)
        #x.append(1) #x[n+1] = 1
        #y.append(1) #y[m+1] = 1
        
        D = np.full((n+1,m+1), np.inf) #D(1:n+1,1:m+1) = inf
        D[0,0] = 0 #D(1,1) = 0
        for i in range(n):
            for j in range(m):
                a1 = D[i,j] + np.power(x[i]-y[j], 2)
                if (i > 1) and (j > 1):
                    if (x[i] > 0) and (y[j] > 0):
                        a1 = D[i,j] + np.power(x[i]-y[j], 2)
                    elif (x[i] < 0) and (y[j] < 0):
                        a1 = D[i,j]
                    elif (x[i] > 0) and (y[j] < 0):
                        a1 = D[i,j] + np.power(x[i], 2) * (-y[j]) #upperbound
                       # a1 = D[i,j] + x[i]^2 #lowerbound
                    elif (x[i] < 0) and (y[j] > 0):
                        a1 = D[i,j] + np.power(y[j], 2) * (-x[i]) #upperbound
                       # a1 = D[i,j] + y[j]^2   #lowerbound
                
                a2 = np.inf
                if (x[i] > 0) and (y[j] > 0):
                    a2 = D[i+1,j] + np.power(x[i]-y[j],2)
                elif (x[i] < 0) and (y[j] < 0):
                    a2 = D[i+1,j]
                elif (x[i] < 0) and (y[j] > 0):
                    a2 = D[i+1,j] + np.power(y[j],2)
                elif (x[i] > 0) and (y[j] < 0):
                    a2 = D[i+1,j] + np.power(x[i],2) * (-y[j]);
                
                a3 = np.inf
                if (x[i] > 0) and (y[j] > 0 ):
                    a3 = D[i,j+1] + np.power(x[i]-y[j],2)
                elif (x[i] < 0) and (y[j] < 0):
                    a3 = D[i,j+1]
                elif (x[i] > 0) and (y[j] < 0):
                    a3 = D[i,j+1] + np.power(x[i],2)
                elif (x[i] < 0) and (y[j] > 0):
                    a3 = D[i,j+1] + np.power(y[j],2) * (-x[i])
                
                D[i+1,j+1] = min([a1, a2, a3])
                
        d = np.sqrt(D[n,m])
        return d#, D

    def _sparse_period(self, start, delta, series, discr = 'minutes'):
        import datetime as dt
        stamps_in_window = []
        for timestamp in series:
            if timestamp >= start:
                if timestamp <= (start+delta):
                    stamps_in_window.append(timestamp)
                else:
                    break
        
        #To sparsed format
        sparsed_stamps = [-int(delta.total_seconds()/60)]
        if stamps_in_window != []:
            sparsed_stamps = [-int((stamps_in_window[0] - start).total_seconds()/60 - 1), 1]
            for i in range(1,len(stamps_in_window)):
                if (stamps_in_window[i] - stamps_in_window[i-1]).total_seconds()/60 > 1:
                        sparsed_stamps.append(-int((stamps_in_window[i] - 
                                                    stamps_in_window[i-1]).total_seconds()/60) - 1)
                sparsed_stamps.append(1)
                
            if start + delta > stamps_in_window[-1]:
                sparsed_stamps.append(-int((start + delta - stamps_in_window[-1]).total_seconds()/60))
            if sparsed_stamps[-1] == 0:
                sparsed_stamps.pop()
        #return stamps_in_window, sparsed_stamps
        return sparsed_stamps
    
    def _positive_area_in_period(self, start, delta, series):
        import datetime as dt
        threshold = dt.timedelta(0, 3*60*60)
        stamps_in_window = []
        for i in range(len(series)):
            if series[i] >= start:
                if series[i] <= (start+delta):
                    stamps_in_window.append(i)
                else:
                    break
        
        positive_time = dt.timedelta(0, 0)
        if stamps_in_window != []:
            start_i = stamps_in_window[0]
            end_i = stamps_in_window[-1]
            if start_i != 0:
                if (series[start_i] - series[start_i-1]) < threshold:
                    positive_time += series[start_i] - start
            if end_i != len(series)-1:
                if series[end_i+1] - series[end_i] < threshold:
                    positive_time += start + delta - series[end_i]
            for i in range(1,end_i - start_i):
                if series[start_i + i] - series[start_i + i - 1] < threshold:
                    positive_time += series[start_i + i] - series[start_i + i - 1]
        
        #return stamps_in_window, sparsed_stamps
        return positive_time / delta

    def single_sliding_window(self, series, kernel, window_start, window_dur, period = 7):
        """
        series - массив (серия) datetime-ов
        period - период повторения ядра в днях
        kernel - ядро в формате спарс-массива. Пример: [-59,1,1,-20,1,-52,...]
        window_start - начало анализируемого окна (timedelta-смещение от 00.00 первого дня серии)
        window_dur - длительность анализируемого окна (в часах, int или float)
        """
        import datetime as dt
        import numpy as np
        #period
        period = dt.timedelta(period)
        window_dur = dt.timedelta(0, window_dur*60*60)
        step = 5 #шаг внутри окна по минутам
        
        initial_point = dt.datetime(series[0].year, series[0].month, series[0].day)
        total_end_point = dt.datetime(series[-1].year, series[-1].month,series[-1].day) + dt.timedelta(1)
        kernel_duration = sum([abs(entry) for entry in kernel])
        
        #(потенциально?) максимальное значение для сравнения
        inactivity_awarp = self._AWarp([-kernel_duration], kernel) 
        if len(kernel) == 1 and sum(kernel)<0:
            #значение на случай поиска "дыр"
            inactivity_awarp = self._AWarp([1]*kernel_duration, kernel) 
        kernel_duration = dt.timedelta(0, kernel_duration*60)
            
        num_of_entries = int((series[-1] - (initial_point + window_start))/ period + 1)
        step = 5 #шаг внутри окна по минутам
        local_distance = []
        for i in range(num_of_entries):
            #начало и конец одного "сезонного" окна
            #внутри должны смещаться ядро и минимизироваться дистанция AWarp
            start_point = i*period + window_start + initial_point
            end_point = i*period + window_start + initial_point + window_dur
                    
            #ход кернеля по окну. сперва в спарс на кажд., потом оптимизировать
            period_from_series = self._sparse_period(start_point, kernel_duration, series)
            min_awarp = self._AWarp(period_from_series, kernel)
            for j in range(1,int((end_point - start_point - kernel_duration).total_seconds()/60/step)):
                delta = dt.timedelta(0,j*60*step)
                cur_awarp = self._AWarp(self._sparse_period(start_point + delta, kernel_duration, series), 
                                        kernel)
                if cur_awarp < min_awarp:
                    min_awarp = cur_awarp
            local_distance.append(min_awarp)
    
        quality_measure = sum((inactivity_awarp - np.array(local_distance))/
                              inactivity_awarp - 0.5)/num_of_entries
        
        return quality_measure
    
    def sliding_area_percentage(self, series, window_start, window_dur, period = 7, negative = False):
        '''
        series - массив (серия) datetime-ов
        period - период повторения ядра в днях
        window_start - начало анализируемого окна (timedelta-смещение от 00.00 первого дня серии)
        window_dur - длительность анализируемого окна (в часах, int или float)
        '''
        import datetime as dt
        import numpy as np
        #period
        period = dt.timedelta(period)
        window_dur = dt.timedelta(0, window_dur*60*60) #вместо шифта и только в одном направлении
        #quality_measures, interperiodic_shifts = [], [] #"качество" по смещениям внутри периода
        step = 5 #шаг внутри окна по минутам
        
        initial_point = dt.datetime(series[0].year, series[0].month, series[0].day)
        total_end_point = dt.datetime(series[-1].year, series[-1].month, series[-1].day) + dt.timedelta(1)
            
        num_of_entries = int((series[-1] - (initial_point + window_start)) / period + 1)
        step = 5 #шаг внутри окна по минутам
        positive_area = []
        for i in range(num_of_entries):
            #начало и конец одного "сезонного" окна. внутри должны смещаться ядро и минимизироваться дистанция AWarp
            start_point = i*period + window_start + initial_point
            end_point = i*period + window_start + initial_point + window_dur
                    
            #ход кернеля по окну. сперва в спарс на кажд., потом оптимизировать
            positive_area.append(self._positive_area_in_period(start_point, window_dur, series))
            
            #print(start_point, end_point, min_awarp)
    
        if negative:
            quality_measure = sum([1-entry for entry in positive_area])/num_of_entries
        else:
            quality_measure = sum(positive_area)/num_of_entries
        #print(quality_measures)
        #if quality_measure == -np.inf:
            #print(num_of_entries, inactivity_awarp)
            #print(local_distance)
        
        return quality_measure
    
    def plot_rose_daily(self, datetimes, margin_type = None, special_day = None):
        """
        margin_type - разделители по часам для разных групп
        special_day - подсветка заданного дня недели для дебага
        """
        import datetime as dt
        import matplotlib.pyplot as plt
        x, y = [], []
        #fig = plt.figure(figsize=(20,10))
        for date in datetimes:
            x.append((date.month - datetimes[0].month)*
                     (dt.date(date.year, date.month + 1, 1) -
                      dt.date(date.year, date.month, 1)).days + 
                      date.day)
            y.append(date.hour + float(date.minute/60))
            if date.weekday() not in [5, 6]:
                plt.plot(x[-1], y[-1], 'r.')
            else:
                plt.plot(x[-1], y[-1], 'b.')
            if special_day != None:
                if date.weekday() == special_day:
                    plt.plot([x[-1],x[-1]], [9,21],'g')
        if margin_type == 'regular':
            plt.plot([x[0],x[-1]], [6.5, 6.5], 'k')
            plt.plot([x[0],x[-1]], [12.5, 12.5], 'k')
            plt.plot([x[0],x[-1]], [17, 17], 'k')
            plt.plot([x[0],x[-1]], [21, 21], 'k')
        elif margin_type == 'night':
            plt.plot([x[0],x[-1]], [9, 9], 'k')
            plt.plot([x[0],x[-1]], [21, 21], 'k')
    
        plt.show()

class TSKernelAPI():
    def __init__(self):
        pass
    
    def check_as_regular(self, series, plot = False):
        import datetime as dt
        methods = SlidingAwarp()
        morning_awarp, evening_awarp, day_awarp, night_awarp = 0, 0, 0, 0
        morning_area, evening_area, day_area, night_area = 0, 0, 0, 0
        for i in range(7):
            morning_awarp += methods.single_sliding_window(series, [-10,1,-5,1,-5,1,-10], dt.timedelta(i, 6.5*60*60), 6)
            evening_awarp += methods.single_sliding_window(series, [-10,1,-5,1,-5,1,-10], dt.timedelta(i, 17*60*60), 4)
            day_awarp += single_sliding_window(track_datetimes, [-260], dt.timedelta(i1, 12.5*60*60), 4.5)
            night_awarp += single_sliding_window(track_datetimes, [-260], dt.timedelta(i1, 22*60*60), 8.5)
            
            morning_area += methods.sliding_area_percentage(series, dt.timedelta(i, 6.5*60*60), 6)
            evening_area += methods.sliding_area_percentage(series, dt.timedelta(i, 17*60*60), 4)
            day_area += methods.sliding_area_percentage(series, dt.timedelta(i, 12.5*60*60), 4.5, negative = True)
            night_area += methods.sliding_area_percentage(series, dt.timedelta(i, 22*60*60), 8.5, negative = True)
            
        morning_awarp = morning_awarp/7
        evening_awarp = evening_awarp/7
        day_awarp = day_awarp/7
        night_awarp = night_awarp/7
        
        morning_area = morning_area/7
        evening_area = evening_area/7
        day_area = day_area/7
        night_area = night_area/7
        if plot:
            methods.plot_rose_daily(series, margin_type='regular')
        ret = [[morning_awarp, evening_awarp, day_awarp, night_awarp],
               [morning_area, evening_area, day_area, night_area]]
        return ret
    
    def check_as_nightbird(self, series, plot = False):
        import datetime as dt
        day_awarp, night_awarp, day_area, night_area = 0, 0, 0, 0
        methods = SlidingAwarp()
        for i in range(7):
            day_awarp += methods.single_sliding_window(series, [-(8*60)], dt.timedelta(i, 9*60*60), 12)
            night_awarp += methods.single_sliding_window(series, [-10,1,-10,1,-10], dt.timedelta(i, 22*60*60), 11)
            day_area += methods.sliding_area_percentage(series, dt.timedelta(i, 9*60*60), 12, negative = True)
            night_area += methods.sliding_area_percentage(series, dt.timedelta(i, 22*60*60), 11)
        day_awarp = day_awarp/7
        night_awarp = night_awarp/7
        day_area = day_area/7
        night_area = night_area/7
        if plot:
            methods.plot_rose_daily(series, margin_type='night')
        return day_awarp, night_awarp, day_area, night_area
    
    