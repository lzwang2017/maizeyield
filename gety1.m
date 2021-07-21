function [ap, y1] = gety1(tpower, uclusters, vmem, RT, y2, y, WA)
  nclusters = length(uclusters);

  x = repmat(y2,1,(1+tpower)*nclusters).*RT;

  nt = 41;
  t01 = [0:1/(nt-1):1]';
  TA = t01.^[0:tpower];
  AP_A = [-TA; TA; -diff(TA); diff(TA)];
  AP_b = [zeros(nt,1); ones(nt,1); 1/40*ones(nt-1,1); 2/40*ones(nt-1,1)];
  AP_L = -inf(1+tpower, 1);
  AP_U = inf(1+tpower, 1);

  ap = zeros(nclusters*(1+tpower), 1);

  for i = 1:nclusters
    fc = find(vmem == uclusters(i));
    if isempty(fc)
      continue
    end
    ap1 = QRGWAb(x(fc, i:nclusters:end), y(fc), WA(fc), AP_A, AP_b, AP_L, AP_U);
    % ap1 = ap1 / max(1, sum(ap1));

    % QP_c = -x(fc, i:nclusters:end)'*(WA(fc).*y(fc));
    % QP_Q = x(fc, i:nclusters:end)'*(repmat(WA(fc),1,11).*x(fc, i:nclusters:end));
    if isempty(ap1)
      error('ap1 error')
    end

    ap(i:nclusters:end) = ap1;
%     if i == 14
%         i
%     end
  end

  ap = reshape(ap, nclusters, 1+tpower);
  f0 = find(sum(abs(ap),2) < 1e-4);
  f1 = find(sum(abs(ap),2) >= 1e-4);
  nf0 = length(f0);
  if nf0 >= 1
    ap(f0,:) = repmat(WA(f1)'*ap(f1,:)/sum(WA(f1)), nf0, 1);
  end
  ap = ap(:);

  y1 = RT * ap;

end
