open Core.Std
open Qptypes


type pub_state =
| Waiting
| Running of string
| Stopped

let pub_state_of_string = function
| "Waiting"  -> Waiting
| "Stopped"  -> Stopped
| s -> Running s

let string_of_pub_state = function
| Waiting  -> "Waiting"
| Stopped  -> "Stopped"
| Running s -> s



type t =
{
    queue           : Queuing_system.t ;
    state           : Message.State.t option ;
    address_tcp     : Address.Tcp.t option ; 
    address_inproc  : Address.Inproc.t option ;
    psi             : Message.Psi.t option;
    progress_bar    : Progress_bar.t option ;
    running         : bool;
}



let debug_env =
  match Sys.getenv "QP_TASK_DEBUG" with
  | Some x -> x <> ""
  | None   -> false


let debug str =
  if debug_env then
    Printf.printf "TASK : %s%!" str



let zmq_context =
  ZMQ.Context.create ()


let bind_socket ~socket_type ~socket ~port =
  let rec loop = function
  | 0 -> failwith @@ Printf.sprintf
        "Unable to bind the %s socket to port : %d "
        socket_type port
  | -1 -> ()
  | i -> 
      try
        ZMQ.Socket.bind socket @@ Printf.sprintf "tcp://*:%d" port;
        loop (-1)
      with
      | Unix.Unix_error _ -> (Time.pause @@ Time.Span.of_float 1. ; loop (i-1) )
      | other_exception -> raise other_exception
  in loop 60;
  ZMQ.Socket.bind socket @@ Printf.sprintf "ipc:///tmp/qp_run:%d" port


let hostname = lazy (
    try
      Unix.gethostname ()
    with
    | _ -> "localhost"
)


let ip_address = lazy (
  match Sys.getenv "QP_NIC" with
  | None ->
      begin
        try
          Lazy.force hostname
          |> Unix.Inet_addr.of_string_or_getbyname
          |> Unix.Inet_addr.to_string
        with
        | Unix.Unix_error _ ->
            failwith "Unable to find IP address from host name."
      end
  | Some interface ->
      begin
        try
          ok_exn Linux_ext.get_ipv4_address_for_interface interface
        with
        | Unix.Unix_error _ ->
            Lazy.force hostname
            |> Unix.Inet_addr.of_string_or_getbyname
            |> Unix.Inet_addr.to_string
      end
)


let reply_ok rep_socket =
    Message.Ok_msg.create ()
    |> Message.Ok_msg.to_string
    |> ZMQ.Socket.send rep_socket 

let reply_wrong_state rep_socket =
    Printf.printf "WRONG STATE\n%!";
    Message.Error_msg.create "Wrong state"
    |> Message.Error_msg.to_string
    |> ZMQ.Socket.send rep_socket 



let stop ~port =
    debug "STOP";
    let req_socket =
      ZMQ.Socket.create zmq_context ZMQ.Socket.req
    and address =
      Printf.sprintf "ipc:///tmp/qp_run:%d" port
    in
    ZMQ.Socket.set_linger_period req_socket 1_000_000;
    ZMQ.Socket.connect req_socket address;

    Message.Terminate (Message.Terminate_msg.create ())
    |> Message.to_string
    |> ZMQ.Socket.send req_socket ;

    let msg =
      ZMQ.Socket.recv req_socket
      |> Message.of_string
    in
    let () =
      match msg with
      | Message.Ok _ -> ()
      | _ -> failwith "Problem in termination"
    in
    ZMQ.Socket.set_linger_period req_socket 1_000;
    ZMQ.Socket.close req_socket


let new_job msg program_state rep_socket pair_socket =

    let state = 
        msg.Message.Newjob_msg.state
    in

    let progress_bar =
        Progress_bar.init
          ~start_value:0.
          ~end_value:1.
          ~bar_length:20
          ~title:(Message.State.to_string state)
    in

    let result = 
      { program_state with
        state           = Some state ;
        progress_bar    = Some progress_bar ;
        address_tcp     = Some msg.Message.Newjob_msg.address_tcp;
        address_inproc  = Some msg.Message.Newjob_msg.address_inproc;
      }
    in
    reply_ok rep_socket;
    string_of_pub_state Waiting
    |> ZMQ.Socket.send pair_socket ;
    result

let change_pub_state msg program_state rep_socket pair_socket =
  let msg = 
    match msg with
    | `Waiting -> Waiting 
    | `Stopped -> Stopped
    | `Running ->
      begin
	let state =
	   match program_state.state with
	   | Some x -> x
	   | None  -> failwith "Trying to change pub state while no job is ready"
        in
        Running (Message.State.to_string state)
      end
  in
  reply_ok rep_socket;
  string_of_pub_state msg
  |> ZMQ.Socket.send pair_socket ;

  program_state

let end_job msg program_state rep_socket pair_socket =

    let failure () =
        reply_wrong_state rep_socket;
        program_state

    and success state =
        reply_ok rep_socket;
        { program_state with
          state        = None ;
          progress_bar = None ;
        }

    in
    match program_state.state with
    | None -> failure ()
    | Some state -> 
      begin
        if (msg.Message.Endjob_msg.state = state) then
          begin
            string_of_pub_state Waiting
            |> ZMQ.Socket.send pair_socket ;
            success state
          end
        else
          failure ()
      end


let connect msg program_state rep_socket =

    let state =
        match program_state.state with
        | Some state -> state
        | None -> assert false
    in
  
    let push_address =
        match msg with
        | Message.Connect_msg.Tcp    -> 
          begin
              match program_state.address_tcp  with
              | Some address -> Address.Tcp address
              | None -> failwith "Error: No TCP address"
          end
        | Message.Connect_msg.Inproc -> 
          begin
              match program_state.address_inproc with
              | Some address -> Address.Inproc address
              | None -> failwith "Error: No inproc address"
          end
        | Message.Connect_msg.Ipc    -> assert false
    in
    
    let new_queue, client_id =
        Queuing_system.add_client program_state.queue
    in
    Message.ConnectReply (Message.ConnectReply_msg.create
        ~state:state ~client_id ~push_address)
    |> Message.to_string
    |> ZMQ.Socket.send rep_socket ;
    { program_state with
      queue = new_queue 
    }


let disconnect msg program_state rep_socket =

    let state, client_id =
      msg.Message.Disconnect_msg.state, 
      msg.Message.Disconnect_msg.client_id
    in

    let failure () =
        reply_wrong_state rep_socket;
        program_state
    
    and success () = 

        let new_program_state = 
            { program_state with
              queue = Queuing_system.del_client ~client_id program_state.queue
            }
        in
        Message.DisconnectReply (Message.DisconnectReply_msg.create ~state)
        |> Message.to_string
        |> ZMQ.Socket.send rep_socket ;
        new_program_state

    in

    match program_state.state with
    | None -> assert false
    | Some state' -> 
      begin
        if (state = state') then
          success ()
        else
          failure ()
      end
    
let del_task msg program_state rep_socket =

    let state, task_id =
      msg.Message.DelTask_msg.state, 
      msg.Message.DelTask_msg.task_id
    in

    let failure () =
        reply_wrong_state rep_socket;
        program_state
    
    and success () = 

        let new_program_state = 
            { program_state with
              queue = Queuing_system.del_task ~task_id program_state.queue
            }
        in
        let more = 
            (Queuing_system.number_of_tasks new_program_state.queue > 0)
        in
        Message.DelTaskReply (Message.DelTaskReply_msg.create ~task_id ~more)
        |> Message.to_string
        |> ZMQ.Socket.send ~block:true rep_socket ; (** /!\ Has to be blocking *)
        new_program_state

    in

    match program_state.state with
    | None -> assert false
    | Some state' -> 
      begin
        if (state = state') then
          success ()
        else
          failure ()
      end



let add_task msg program_state rep_socket =

    let state, task =
        msg.Message.AddTask_msg.state,
        msg.Message.AddTask_msg.task
    in

    let increment_progress_bar = function
      | Some bar -> Some (Progress_bar.increment_end bar)
      | None -> None
    in

    let rec add_task_triangle program_state imax = function
      | 0 -> program_state
      | i -> 
          let task = 
             Printf.sprintf "%d %d" i imax
          in
          let new_program_state =
             { program_state with
               queue = Queuing_system.add_task ~task program_state.queue ;
               progress_bar = increment_progress_bar program_state.progress_bar ;
             }
          in
          add_task_triangle new_program_state imax (i-1)
    in
    
    let rec add_task_range program_state i = function
      | j when (j < i) -> program_state
      | j -> 
          let task = 
             Printf.sprintf "%d" j
          in
          let new_program_state =
             { program_state with
               queue = Queuing_system.add_task ~task program_state.queue ;
               progress_bar = increment_progress_bar program_state.progress_bar ;
             }
          in
          add_task_range new_program_state i (j-1)
    in
    
    let new_program_state = function
      | "triangle" :: i_str :: [] ->
          let imax = 
              Int.of_string i_str
          in
          add_task_triangle program_state imax imax
      | "range" :: i_str :: j_str :: [] ->
          let i, j = 
              Int.of_string i_str,
              Int.of_string j_str
          in
          add_task_range program_state i j
      | _ -> 
          { program_state with
            queue = Queuing_system.add_task ~task program_state.queue ;
            progress_bar = increment_progress_bar program_state.progress_bar ;
          }
    in 

    let result = 
        String.split ~on:' ' task
        |> List.filter ~f:(fun x -> x <> "")
        |> new_program_state 
    in
    reply_ok rep_socket;
    result



let get_task msg program_state rep_socket pair_socket =

    let state, client_id =
      msg.Message.GetTask_msg.state, 
      msg.Message.GetTask_msg.client_id
    in

    let failure () =
        reply_wrong_state rep_socket;
        program_state
    
    and success () = 

        let new_queue, task_id, task =
            Queuing_system.pop_task ~client_id program_state.queue
        in
        if (Queuing_system.number_of_queued new_queue = 0) then
          string_of_pub_state Waiting 
          |> ZMQ.Socket.send pair_socket
        else
          string_of_pub_state (Running (Message.State.to_string state))
          |> ZMQ.Socket.send pair_socket;

        let new_program_state = 
            { program_state with
              queue = new_queue
            }
        in

        Message.GetTaskReply (Message.GetTaskReply_msg.create ~task ~task_id)
        |> Message.to_string
        |> ZMQ.Socket.send rep_socket ;
        new_program_state

    in

    match program_state.state with
    | None -> assert false
    | Some state' -> 
      begin
        if (state = state') then
          success ()
        else
          failure ()
      end



let task_done msg program_state rep_socket =

    let state, client_id, task_id =
        msg.Message.TaskDone_msg.state,
        msg.Message.TaskDone_msg.client_id,
        msg.Message.TaskDone_msg.task_id
    in

    let increment_progress_bar = function
      | Some bar -> Some (Progress_bar.increment_cur bar)
      | None -> None
    in

    let failure () =
        reply_wrong_state rep_socket;
        program_state
    
    and success () = 
        let result = 
          { program_state with
            queue = Queuing_system.end_task ~task_id ~client_id program_state.queue ;
            progress_bar = increment_progress_bar program_state.progress_bar ;
          }
        in
        reply_ok rep_socket;
        result
    in

    match program_state.state with
    | None -> assert false
    | Some state' -> 
      begin
        if (state = state') then
          success ()
        else
          failure ()
      end


let put_psi msg rest_of_msg program_state rep_socket =

    let psi_local =
        match msg.Message.PutPsi_msg.psi with
        | Some x -> x
        | None ->
          begin
            let psi_det, psi_coef, energy =
              match rest_of_msg with
              | [ x ; y ; e ] -> x, y, e
              | _ -> failwith "Badly formed put_psi message"
            in
            Message.Psi.create
              ~n_state:msg.Message.PutPsi_msg.n_state
              ~n_det:msg.Message.PutPsi_msg.n_det
              ~psi_det_size:msg.Message.PutPsi_msg.psi_det_size
              ~n_det_generators:msg.Message.PutPsi_msg.n_det_generators
              ~n_det_selectors:msg.Message.PutPsi_msg.n_det_selectors
              ~psi_det
              ~psi_coef
              ~energy
          end
    in
    let new_program_state =
        { program_state with
          psi = Some psi_local
        }
    and client_id =
      msg.Message.PutPsi_msg.client_id
    in
    Message.PutPsiReply (Message.PutPsiReply_msg.create ~client_id)
    |> Message.to_string
    |> ZMQ.Socket.send rep_socket;

    new_program_state


let get_psi msg program_state rep_socket =

      let client_id =
          msg.Message.GetPsi_msg.client_id
      in
      match program_state.psi with
      | None -> failwith "No wave function saved in TaskServer"
      | Some psi -> 
          Message.GetPsiReply (Message.GetPsiReply_msg.create ~client_id ~psi)
          |> Message.to_string_list 
          |> ZMQ.Socket.send_all rep_socket;
      program_state



let terminate program_state rep_socket =
    reply_ok rep_socket;
    { program_state with
      psi = None;
      address_tcp = None;
      address_inproc = None;
      running = false
    }

  
let error msg program_state rep_socket =
    Message.Error (Message.Error_msg.create msg)
    |> Message.to_string
    |> ZMQ.Socket.send rep_socket ;
    program_state

let start_pub_thread ~port = 
    Thread.create (fun () ->
      let timeout =
        1000
      in

      let pair_socket =
        ZMQ.Socket.create zmq_context ZMQ.Socket.pair
      and address = 
        "inproc://pair" 
      in
      ZMQ.Socket.connect pair_socket address;

      let pub_socket =
        ZMQ.Socket.create zmq_context ZMQ.Socket.pub
      in
      bind_socket ~socket_type:"PUB" ~socket:pub_socket ~port;

      let pollitem =
        ZMQ.Poll.mask_of 
          [| (pair_socket, ZMQ.Poll.In) |]
      in

      let rec run state = 
        let new_state =
          let polling = 
            ZMQ.Poll.poll ~timeout pollitem
          in
          if (polling.(0) = Some ZMQ.Poll.In) then
            ZMQ.Socket.recv ~block:false pair_socket
            |> pub_state_of_string
          else
            state
        in
        ZMQ.Socket.send pub_socket @@ string_of_pub_state new_state;
        match state with
        | Stopped -> ()
        | _ -> run new_state
      in
      run Waiting;
      ZMQ.Socket.set_linger_period pair_socket 1000 ;
      ZMQ.Socket.close pair_socket;
      ZMQ.Socket.set_linger_period pub_socket 1000 ;
      ZMQ.Socket.close pub_socket;
    )

let run ~port =

    (** Bind inproc socket for changing state of pub *)
    let pair_socket =
      ZMQ.Socket.create zmq_context ZMQ.Socket.pair
    and address = 
      "inproc://pair" 
    in
    ZMQ.Socket.bind pair_socket address;

    let pub_thread = 
      start_pub_thread ~port:(port+1) ()
    in

    (** Bind REP socket *)
    let rep_socket =
      ZMQ.Socket.create zmq_context ZMQ.Socket.rep
    in
    ZMQ.Socket.set_linger_period rep_socket 1_000_000;
    bind_socket "REP" rep_socket port;

    let initial_program_state =
    {   queue = Queuing_system.create () ;
        running = true ;
        psi = None;
        state = None;
        address_tcp = None;
        address_inproc = None;
        progress_bar = None ;
    }
    in

    (** ZMR polling item *)
    let pollitem =
      ZMQ.Poll.mask_of
        [| (rep_socket, ZMQ.Poll.In) |]
    in

    let address =
      Printf.sprintf "tcp://%s:%d" (Lazy.force ip_address) port
    in
    Printf.printf "Task server running : %s\n%!" address;


    (** Main loop *)
    let rec main_loop program_state = function
    | false -> ()
    | true ->
        let polling =
            ZMQ.Poll.poll ~timeout:1000 pollitem
        in
        if (polling.(0) <> Some ZMQ.Poll.In) then
          main_loop program_state true
        else
          begin
              let program_state = 
                  match program_state.progress_bar with
                  | None -> program_state
                  | Some bar -> 
                    if bar.Progress_bar.dirty then 
                        { program_state with
                          progress_bar = Some (Progress_bar.display bar)
                        }
                    else
                        program_state
              in

              (** Extract message *)
              let raw_message, rest =
                  match ZMQ.Socket.recv_all rep_socket with
                  | x :: rest -> x, rest
                  | [] -> failwith "Badly formed message"
              in
              let message =
                  Message.of_string raw_message
              in

              (** Debug input *)
              Printf.sprintf "q:%d  r:%d  n:%d  : %s\n%!"
              (Queuing_system.number_of_queued  program_state.queue)
              (Queuing_system.number_of_running program_state.queue)
              (Queuing_system.number_of_tasks   program_state.queue)
              (Message.to_string message)
              |> debug;

              let new_program_state = 
                try
                  match program_state.state, message with
                  | _     , Message.Terminate   _ -> terminate program_state rep_socket
                  | _     , Message.PutPsi      x -> put_psi x rest program_state rep_socket
                  | _     , Message.GetPsi      x -> get_psi x program_state rep_socket
                  | None  , Message.Newjob      x -> new_job x program_state rep_socket pair_socket
                  | _     , Message.Newjob      _ -> error "A job is already running" program_state rep_socket
                  | Some _, Message.Endjob      x -> end_job x program_state rep_socket pair_socket
                  | Some _, Message.SetRunning    -> change_pub_state `Running program_state rep_socket pair_socket
                  | _, Message.SetWaiting    -> change_pub_state `Waiting program_state rep_socket pair_socket
                  | _, Message.SetStopped    -> change_pub_state `Stopped program_state rep_socket pair_socket
                  | None  ,                     _ -> error "No job is running" program_state rep_socket
                  | Some _, Message.Connect     x -> connect x program_state rep_socket
                  | Some _, Message.Disconnect  x -> disconnect x program_state rep_socket
                  | Some _, Message.AddTask     x -> add_task x program_state rep_socket
                  | Some _, Message.DelTask     x -> del_task x program_state rep_socket
                  | Some _, Message.GetTask     x -> get_task x program_state rep_socket pair_socket
                  | Some _, Message.TaskDone    x -> task_done x program_state rep_socket
                  | _     , _                     ->
                    error ("Invalid message : "^(Message.to_string message))  program_state rep_socket
                with
                | Failure f -> 
                    error (f^" : "^raw_message) program_state rep_socket
                | Assert_failure (f,i,j) ->
                    error (Printf.sprintf "%s:%d:%d : %s" f i j raw_message) program_state rep_socket

              in
              main_loop new_program_state new_program_state.running
          end
    in main_loop initial_program_state true;

    ZMQ.Socket.send pair_socket @@ string_of_pub_state Stopped;
    Thread.join pub_thread;
    ZMQ.Socket.close rep_socket




    
