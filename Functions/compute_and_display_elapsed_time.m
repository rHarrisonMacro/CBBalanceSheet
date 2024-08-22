function compute_and_display_elapsed_time(startTime,finishTime)

processingTime = finishTime - startTime;
days = processingTime(3);
hours = processingTime(4);
mins = processingTime(5);
secs = processingTime(6);
if secs<0
    secs = 60 + secs;
    mins = mins - 1;
end
if mins<0
    mins = mins + 60;
    hours = hours - 1;
end
if hours<0
    hours = hours + 24;
    days = days - 1;
end
disp(['Total computing time = ' num2str(days) ' days, ' ...
    num2str(hours) ' hours, ' num2str(mins) ' minutes, ' ...
    num2str(secs) ' seconds.']);

end

